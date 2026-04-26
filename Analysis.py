import pandas as pd
import glob
import os

# Get all CSV files inside outputs directory
files = glob.glob(os.path.join("outputs", "*_STO-*G.csv"))

results = []

for file in files:
      
    # Extract molecule and basis set from filename
    filename = os.path.basename(file)
    name_parts = filename.replace('.csv', '').split('_')
    molecule = name_parts[0]
    basis_set = name_parts[1]
    
    # Read the CSV
    df = pd.read_csv(file)
    
    # Identify columns
    bond_col = df.columns[0]
    energy_col = df.columns[1]
    
    # 1. Equilibrium bond length (r_e) and Minimum Energy (E_min)
    min_idx = df[energy_col].idxmin()
    e_min = df.loc[min_idx, energy_col]
    r_e = df.loc[min_idx, bond_col]
    
    # 2. Dissociation Energy Limit (E_dissociation)
    # Approximated by the energy at the maximum calculated distance in the dataset
    max_bond_idx = df[bond_col].idxmax()
    e_diss = df.loc[max_bond_idx, energy_col]
    
    # 3. Binding Energy (D_e) = E_dissociation - E_min
    de = e_diss - e_min
    
    results.append({
        'Molecule': molecule,
        'Basis_Set': basis_set,
        'r_e (Angstrom)': r_e,
        'E_min (Hartree)': e_min,
        'E_dissociation (Hartree)': e_diss,
        'D_e (Hartree)': de
    })

# Convert to DataFrame and sort
results_df = pd.DataFrame(results).sort_values(by=['Molecule', 'Basis_Set']).reset_index(drop=True)

# Save to CSV
results_df.to_csv('extracted_properties.csv', index=False)
print(results_df.to_markdown())