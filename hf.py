"""
Restricted Hartree-Fock (RHF) — Bond Dissociation Curve

Computes the RHF total energy of H2 or HeH+ across a range of
bond lengths and writes the results to a CSV for plotting.
Molecules: H2 or HeH+
Basis sets: STO-2G to STO-6G
"""

import math
import csv
import logging
import numpy as np
from datetime import datetime

# HARDCODED VALUES
DEFAULT_MOLECULE   = "H2"                   # "H2" or "HeH+"
DEFAULT_BASIS_N    = 3                      # n in STO-nG, integer 1-6
DEFAULT_R_START    = 0.3                    # bond length scan start, Angstrom
DEFAULT_R_END      = 4.0                    # bond length scan end, Angstrom
DEFAULT_STEP_SIZE  = 0.05                   # step size between bond lengths, Angstrom
DEFAULT_OUTPUT     = "dissociation.csv"     #Output file name
ENABLE_LOGGING     = True                   # log file generation

# CONSTANTS
BOHR = 1.8897259886   # 1 Angstrom in Bohr (atomic units of length)
NUCLEAR_CHARGE = {"H": 1, "He": 2}

# STO-nG BASIS SET DATA
# Each entry is a list of (exponent, contraction_coefficient) pairs for the
# 1s Slater-type orbital expansion.  
# Values from the https://www.basissetexchange.org/
# The STO-nG basis approximates a single Slater 1s orbital as a sum of n
# Gaussian functions:  chi(r) = sum_k c_k * g_k(r; alpha_k)

BASIS_DATA = {
    "H": {
        2: [(1.309756377, 0.4301284983), (0.2331359749, 0.6789135305)],
        3: [(3.425250914, 0.1543289673), (0.6239137298, 0.5353281423), (0.1688554040, 0.4446345422)],
        4: [(8.021420155, 0.05675242080), (1.467821061, 0.2601413550), (0.4077767635, 0.5328461143), (0.1353374420, 0.2916254405)],
        5: [(17.38354739, 0.02214055312), (3.185489246, 0.1135411520), (0.8897299079, 0.3318161484), (0.3037874103, 0.4825700713), (0.1144784984, 0.1935721966)],
        6: [(35.52322122, 0.009163596281), (6.513143725, 0.04936149294), (1.822142904, 0.1685383049), (0.6259552659, 0.3705627997), (0.2430767471, 0.4164915298), (0.1001124280, 0.1303340841)],
    },
    "He": {
        2: [(2.432879285, 0.4301284983), (0.4330512863, 0.6789135305)],
        3: [(6.362421394, 0.1543289673), (1.158922999, 0.5353281423), (0.3136497915, 0.4446345422)],
        4: [(14.89982967, 0.05675242080), (2.726485258, 0.2601413550), (0.7574474599, 0.5328461143), (0.2513900027, 0.2916254405)],
        5: [(32.29002972, 0.02214055312), (5.917062849, 0.1135411520), (1.652677933, 0.3318161484), (0.5642866953, 0.4825700713), (0.2126444063, 0.1935721966)],
        6: [(65.98456824, 0.009163596281), (12.09819836, 0.04936149294), (3.384639924, 0.1685383049), (1.162715163, 0.3705627997), (0.4515163224, 0.4164915298), (0.1859593559, 0.1303340841)],
    },
}

# BASIS SET CONSTRUCTION

def normalization(alpha):
    # Normalisation for a 1s Gaussian so that <g|g> = 1.
    # N(alpha) = (2*alpha/pi)^(3/4)
    return (2.0 * alpha / math.pi) ** 0.75


def build_basis(elem_A, pos_A, elem_B, pos_B, n):
    # One contracted basis function per atom.
    # Each function is a list of primitive Gaussian dicts.
    # Input positions are in Angstrom; we convert to Bohr here.
    basis = []
    for elem, pos_ang in [(elem_A, pos_A), (elem_B, pos_B)]:
        center = np.array(pos_ang) * BOHR        # convert to atomic units
        primitives = []
        for alpha, coeff in BASIS_DATA[elem][n]:
            primitives.append({
                "alpha":  alpha,
                "coeff":  coeff,
                "norm":   normalization(alpha),
                "center": center,
            })
        basis.append(primitives)
    return basis

# ONE AND TWO ELECTRON INTEGRALS

def boys0(x):
    # F_0(x) = (1/2) sqrt(pi/x) erf(sqrt(x)), limiting value 1 as x->0.
    # This comes up in every Coulomb-type integral (nuclear attraction, ERI).
    if x < 1e-8:
        return 1.0 - x / 3.0
    return 0.5 * math.sqrt(math.pi / x) * math.erf(math.sqrt(x))


def gaussian_product_center(ai, ri, aj, rj):
    # The product of two Gaussians is a Gaussian at the weighted centre rp.
    # Returns the combined exponent p and the new centre rp.
    p  = ai + aj
    rp = (ai * ri + aj * rj) / p
    return p, rp


def overlap_primitive(gi, gj):
    # <g_i|g_j> for two 1s Gaussian primitives.
    p = gi["alpha"] + gj["alpha"]
    q = gi["alpha"] * gj["alpha"] / p
    r2 = np.dot(gi["center"] - gj["center"], gi["center"] - gj["center"])
    return gi["norm"] * gj["norm"] * (math.pi / p)**1.5 * math.exp(-q * r2)


def kinetic_primitive(gi, gj):
    # <g_i | -1/2 nabla^2 | g_j> for two 1s Gaussian primitives.
    # Simplifies to q*(3 - 2*q*r^2) times the overlap, where q = ai*aj/(ai+aj).
    p  = gi["alpha"] + gj["alpha"]
    q  = gi["alpha"] * gj["alpha"] / p
    r2 = np.dot(gi["center"] - gj["center"], gi["center"] - gj["center"])
    S  = gi["norm"] * gj["norm"] * (math.pi / p)**1.5 * math.exp(-q * r2)
    return q * (3.0 - 2.0 * q * r2) * S


def nuclear_primitive(gi, gj, rc, Zc):
    # <g_i | -Z/|r-Rc| | g_j>.  The Gaussian product theorem collapses the
    # two-centre product to one Gaussian at rp, and the 1/r part becomes Boys F_0.
    p, rp = gaussian_product_center(gi["alpha"], gi["center"], gj["alpha"], gj["center"])
    q  = gi["alpha"] * gj["alpha"] / p
    r2 = np.dot(gi["center"] - gj["center"], gi["center"] - gj["center"])
    rpc2 = np.dot(rp - rc, rp - rc)
    prefactor = gi["norm"] * gj["norm"] * (-Zc) * (2.0 * math.pi / p) * math.exp(-q * r2)
    return prefactor * boys0(p * rpc2)


def eri_primitive(gi, gj, gk, gl):
    # (ij|kl) two-electron repulsion integral.
    # Both pairs collapse via the Gaussian product theorem, 1/r12 becomes Boys F_0.
    pij, rij = gaussian_product_center(gi["alpha"], gi["center"], gj["alpha"], gj["center"])
    pkl, rkl = gaussian_product_center(gk["alpha"], gk["center"], gl["alpha"], gl["center"])

    q_ij = gi["alpha"] * gj["alpha"] / pij
    q_kl = gk["alpha"] * gl["alpha"] / pkl
    r2_ij = np.dot(gi["center"] - gj["center"], gi["center"] - gj["center"])
    r2_kl = np.dot(gk["center"] - gl["center"], gk["center"] - gl["center"])
    rpq2  = np.dot(rij - rkl, rij - rkl)

    # combined exponent for the Boys function argument
    theta = pij * pkl / (pij + pkl)

    prefactor = (
        2.0 * math.pi**2 / (pij * pkl)
        * math.sqrt(math.pi / (pij + pkl))
        * math.exp(-q_ij * r2_ij - q_kl * r2_kl)
    )

    norm_product  = gi["norm"]  * gj["norm"]  * gk["norm"]  * gl["norm"]
    coeff_product = gi["coeff"] * gj["coeff"] * gk["coeff"] * gl["coeff"]

    return norm_product * coeff_product * prefactor * boys0(theta * rpq2)


def contracted_integral(func, ci, cj, *args):
    # Sum the given primitive integral over all primitive pairs in ci and cj.
    total = 0.0
    for gi in ci:
        for gj in cj:
            total += gi["coeff"] * gj["coeff"] * func(gi, gj, *args)
    return total


def build_overlap_matrix(basis):
    # S[i,j] = <phi_i | phi_j>
    n = len(basis)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = contracted_integral(overlap_primitive, basis[i], basis[j])
            S[i, j] = val
            S[j, i] = val   # S is symmetric
    return S


def build_kinetic_matrix(basis):
    # T[i,j] = <phi_i | -1/2 nabla^2 | phi_j>
    n = len(basis)
    T = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = contracted_integral(kinetic_primitive, basis[i], basis[j])
            T[i, j] = val
            T[j, i] = val
    return T


def build_nuclear_matrix(basis, nuclei):
    # V[i,j] = sum over nuclei of <phi_i | -Zk/|r-Rk| | phi_j>
    n = len(basis)
    V = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            # sum contributions from all nuclei
            val = sum(
                contracted_integral(nuclear_primitive, basis[i], basis[j], rc, Zc)
                for rc, Zc in nuclei
            )
            V[i, j] = val
            V[j, i] = val
    return V


def build_eri_tensor(basis):
    # 4-index ERI tensor.  8-fold symmetry means (ij|kl) = (ji|kl) = (ij|lk) = (kl|ij)...
    # so we only compute unique index combinations and fill the rest by symmetry.
    n = len(basis)
    G = np.zeros((n, n, n, n))

    for i in range(n):
        for j in range(i, n):
            for k in range(n):
                for l in range(k, n):
                    val = 0.0
                    # sum over all combinations of primitives in the four functions
                    for gi in basis[i]:
                        for gj in basis[j]:
                            for gk in basis[k]:
                                for gl in basis[l]:
                                    val += eri_primitive(gi, gj, gk, gl)
                    # fill all 8 symmetry-equivalent entries at once
                    G[i,j,k,l] = G[j,i,k,l] = G[i,j,l,k] = G[j,i,l,k] = val
                    G[k,l,i,j] = G[l,k,i,j] = G[k,l,j,i] = G[l,k,j,i] = val
    return G




# SCF 

def build_fock_matrix(H_core, D, G):
    # F = H_core + J - 0.5*K
    # J is the Coulomb term, K is the exchange term.  einsum keeps this concise.
    J = np.einsum("kl,ijkl->ij", D, G)
    K = np.einsum("kl,ilkj->ij", D, G)
    return H_core + J - 0.5 * K


def build_density_matrix(C, n_elec):
    # D = 2 * C_occ @ C_occ^T.  Factor of 2 for spin (alpha + beta).
    # Odd electron count: the last occupied orbital gets weight 1 instead of 2.
    n_occ = n_elec // 2
    D = 2.0 * C[:, :n_occ] @ C[:, :n_occ].T
    if n_elec % 2 == 1:
        # singly occupied orbital — one unpaired electron
        D += np.outer(C[:, n_occ], C[:, n_occ])
    return D


def electronic_energy(D, H_core, F):
    # E_elec = 0.5 * Tr[D * (H_core + F)].  The 0.5 avoids double-counting e-e repulsion.
    return 0.5 * np.einsum("ij,ij->", D, H_core + F)


def nuclear_repulsion_energy(nuclei):
    # E_nuc = sum_{A<B} Z_A*Z_B / R_AB.  Constant for a fixed geometry.
    E = 0.0
    for a in range(len(nuclei)):
        for b in range(a + 1, len(nuclei)):
            ra, Za = nuclei[a]
            rb, Zb = nuclei[b]
            E += Za * Zb / np.linalg.norm(ra - rb)
    return E


def run_scf(H_core, S, G, n_elec, max_iter=100, tol=1e-8):
    # SCF loop: build F -> solve eigenvalue problem -> update D -> check convergence.
    # We work in the orthogonal basis (X = S^{-1/2}) so we can use a plain eigh call.

    # X = S^{-1/2}: transforms AOs to an orthonormal basis
    s_vals, s_vecs = np.linalg.eigh(S)
    X = s_vecs @ np.diag(s_vals ** -0.5) @ s_vecs.T

    # Zero density = core Hamiltonian guess (no electron-electron term on first iteration)
    D = np.zeros((len(H_core), len(H_core)))
    E_prev = 0.0

    for iteration in range(1, max_iter + 1):

        # H atom / one-electron case: Fock matrix is just H_core
        if n_elec == 1:
            F = H_core.copy()
        else:
            F = build_fock_matrix(H_core, D, G)

        # Diagonalise in orthogonal basis, then back-transform coefficients
        F_prime = X.T @ F @ X
        orbital_energies, C_prime = np.linalg.eigh(F_prime)
        C = X @ C_prime

        D_new  = build_density_matrix(C, n_elec)
        E_new  = electronic_energy(D_new, H_core, F)
        delta_E = abs(E_new - E_prev)

        logging.info(f"    iter {iteration:3d}   E_elec = {E_new:16.10f}   dE = {delta_E:.3e}")

        if iteration > 1 and delta_E < tol:
            logging.info("    SCF converged.")
            return E_new

        D, E_prev = D_new, E_new

    logging.warning("    SCF did not converge!")
    return E_prev

# MAIN CALCULATIONS

def compute_energy(molecule, n, bond_length_ang):
    # Put atom A at the origin and atom B along z at the given bond length.
    # Returns E_total = E_elec + E_nuc in Hartree.

    # Define atoms and molecular charge
    if molecule == "H2":
        elem_A, elem_B, charge = "H", "H", 0
    else:  # HeH+
        elem_A, elem_B, charge = "He", "H", 1

    pos_A = [0.0, 0.0, 0.0]
    pos_B = [0.0, 0.0, bond_length_ang]

    n_elec = NUCLEAR_CHARGE[elem_A] + NUCLEAR_CHARGE[elem_B] - charge

    # Nuclear positions in Bohr for integral evaluation
    nuclei = [
        (np.array(pos_A) * BOHR, NUCLEAR_CHARGE[elem_A]),
        (np.array(pos_B) * BOHR, NUCLEAR_CHARGE[elem_B]),
    ]

    # Build all integrals
    basis  = build_basis(elem_A, pos_A, elem_B, pos_B, n)
    S      = build_overlap_matrix(basis)
    T      = build_kinetic_matrix(basis)
    V      = build_nuclear_matrix(basis, nuclei)
    H_core = T + V
    G      = build_eri_tensor(basis)
    E_nuc  = nuclear_repulsion_energy(nuclei)

    E_elec = run_scf(H_core, S, G, n_elec)
    return E_elec + E_nuc




#INPUT
def get_inputs():
    # Collect run settings.  Pressing Enter on any prompt uses the default from the top.

    print()
    print("Select Molecule:  a) H2       b) HeH+")
    mol_input = input("                  > ").strip().lower()
    if mol_input == "b" or mol_input == "HeH+":
        molecule = "HeH+"
    elif mol_input == "a" or mol_input == "H2":
        molecule = "H2"
    else:
        print(f"  Unrecognised, using default: {DEFAULT_MOLECULE}")
        molecule = DEFAULT_MOLECULE

    print()
    print("Basis set STO-nG:   n = ", end="")
    n_input = input().strip()
    try:
        n = int(n_input)
        if n not in range(2, 7):
            raise ValueError
    except ValueError:
        print(f"  Invalid, using default: n={DEFAULT_BASIS_N}")
        n = DEFAULT_BASIS_N

    print()
    print("Bond Length Range")
    print("  Start: ", end="")
    try:
        r_start = float(input().strip())
    except ValueError:
        r_start = DEFAULT_R_START

    print("  End:   ", end="")
    try:
        r_end = float(input().strip())
    except ValueError:
        r_end = DEFAULT_R_END

    print("  Step:  ", end="")
    try:
        step = float(input().strip())
        if step <= 0:
            raise ValueError
    except ValueError:
        step = DEFAULT_STEP_SIZE
        print(f"  Bad value, using default step: {DEFAULT_STEP_SIZE} Ang")

    print()
    print(f"Output [{DEFAULT_OUTPUT}]: ", end="")
    out_input = input().strip()
    out_file = out_input if out_input else DEFAULT_OUTPUT

    return molecule, n, r_start, r_end, step, out_file



def main():
    if ENABLE_LOGGING:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file  = f"log_{timestamp}.log"
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format="%(asctime)s  %(levelname)s  %(message)s",
            datefmt="%H:%M:%S",
        )
    else:
        # Suppress all logging calls (CRITICAL and below)
        logging.disable(logging.CRITICAL)
        log_file = "Disabled"
    # Get run parameters from the user
    molecule, n, r_start, r_end, step, out_file = get_inputs()

    basis_label = f"STO-{n}G"
    logging.info(
        f"Run started — molecule={molecule}, basis={basis_label}, "
        f"r=[{r_start}, {r_end}] Ang, step={step} Ang"
    )

    # Generate bond lengths from start to end
    distances = []
    d = r_start
    while d <= r_end + 1e-9:   
        distances.append(round(d, 10))
        d += step

    # CSV column headers
    col_r = "Bond_length_ang"
    col_E = f"STO-{n}G_Energy"

    rows = []
    print()
    for i, d in enumerate(distances, 1):
        print(f"  [{i:3d}/{len(distances)}]  d = {d:.4f} Ang ...", end=" ", flush=True)
        logging.info(f"  d = {d:.6f} Ang")
        try:
            E = compute_energy(molecule, n, d)
            print(f"E = {E:.8f} Ha")
            logging.info(f"    E_total = {E:.10f} Ha")
            rows.append({col_r: round(d, 6), col_E: round(E, 10)})
        except Exception as err:
            print(f"FAILED: {err}")
            logging.error(f"    FAILED at d={d:.6f}: {err}")
            rows.append({col_r: round(d, 6), col_E: ""})

    # Write the two-column CSV
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[col_r, col_E])
        writer.writeheader()
        writer.writerows(rows)

    print()
    print(f"Done.  Results -> {out_file}")
    print(f"       Log     -> {log_file}")

    logging.info(f"Finished. {len(rows)} points written to {out_file}")


if __name__ == "__main__":
    main()