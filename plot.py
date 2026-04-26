import pandas as pd
import matplotlib.pyplot as plt
import os

# Define the parameter space we scanned
molecules = ["H2", "HeH+"]
basis_sets = [2, 3, 4, 5, 6]

# Create a figure with two side-by-side subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, mol in zip(axes, molecules):
    for n in basis_sets:
        filename = f"outputs\\{mol}_STO-{n}G.csv"
        
        # Check if file exists before trying to plot
        if os.path.exists(filename):
            df = pd.read_csv(filename)
            
            # Grab columns by position (0 is distance, 1 is energy)
            distances = df.iloc[:, 0]
            energies = df.iloc[:, 1]
            
            # Plot the curve
            ax.plot(distances, energies, label=f"STO-{n}G", linewidth=2)
        else:
            print(f"Warning: {filename} not found.")

    # Format the subplot
    ax.set_title(f"Bond Dissociation Curve: {mol}", fontsize=14, fontweight="bold")
    ax.set_xlabel("Bond Length (Å)", fontsize=12)
    ax.set_ylabel("Total Energy (Hartree)", fontsize=12)
    ax.legend(title="Basis Set")
    
    # Add a subtle grid for easier reading
    ax.grid(True, linestyle="--", alpha=0.7)
    
    if mol=='H2':
        ax.set_xlim(0.3,3.0)
    else:
        ax.set_xlim(0.5, 3.0)

# Adjust layout to prevent overlapping text
plt.tight_layout()

# Save a high-res image and then display it
plt.savefig("dissociation_curves.png", dpi=300)
print("\nPlot saved as 'dissociation_curves.png'")
plt.show()