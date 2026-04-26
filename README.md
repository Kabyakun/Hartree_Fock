# Hartree Fock Implementation

The folder contains files to compute, analyse, and visualize bond dissociation curves for simple diatomic molecules (H₂ and HeH⁺) using Restricted Hartree-Fock (RHF) calculations written from scratch.

## Project Structure

The workflow is broken down into four main scripts:

* **`hf.py`**: The core computational engine. It performs RHF calculations by computing one- and two-electron integrals and running the Self-Consistent Field (SCF) loop. It evaluates the total energy for H₂ or HeH⁺ over a range of bond lengths using STO-nG basis sets (where n=2 to 6).
* **`run_hf_scans.ps1`**: A PowerShell batch script that automates `hf.py`. It loops through the defined molecules and basis sets, piping simulated inputs into the Python script to generate a complete set of CSV data files (e.g., `H2_STO-3G.csv`).
* **`Analysis.py`**: A data extraction script that parses the generated CSVs. It calculates key molecular properties including Equilibrium Bond Length ($r_e$), Minimum Energy ($E_{min}$), Dissociation Energy Limit ($E_{dissociation}$), and Binding Energy ($D_e$). Results are saved to `extracted_properties.csv`.
* **`plot.py`**: A visualization script. It reads the CSV files and uses `matplotlib` to generate side-by-side plots of the dissociation curves, comparing the accuracy of the different basis sets. The output is saved as `dissociation_curves.png`.

## Requirements

* Python 3.x
* `numpy` (for matrix operations and integral evaluation in `hf.py`)
* `pandas` (for data manipulation in `Analysis.py` and `plot.py`)
* `matplotlib` (for plotting in `plot.py`)

## Usage Instructions

1.  **Generate Data**: 
    Run the PowerShell script to batch generate all dissociation curve CSV files automatically:
    ```powershell
    .\run_hf_scans.ps1
    ```
    *(Alternatively, run `python hf.py` to calculate a specific curve interactively).*

2.  **Analyze Results**:
    Once the CSV files are generated in the directory, extract the key molecular properties:
    ```bash
    python Analysis.py
    ```
    This will print a markdown table to the console and save the data to `extracted_properties.csv`.

3.  **Visualize Curves**:
    Generate the dissociation energy plots:
    ```bash
    python plot.py
    ```
    This will display the graphs and save them as a high-resolution PNG image.
