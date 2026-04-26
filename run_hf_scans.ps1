# batch runner for HF scans across STO-nG basis sets

# Define our molecules
$molecules = @(
    #@{ Key="a"; Name="H2" },
    @{ Key="b"; Name="HeH+" }
)

# Loop through both molecules
foreach ($mol in $molecules) {

    # Loop through each basis sets 
    for ($n = 2; $n -le 6; $n++) {

        # Create a output file the run
        $outFile = "$($mol.Name)_STO-$($n)G.csv"
        Write-Host " Running $($mol.Name) with STO-$($n)G basis set..." -ForegroundColor Cyan
        
        # Create an array of inputs
        $Inputs = @(
            $mol.Key,      # Select Molecule (a or b)
            $n.ToString(), # Basis set (2-6)
            "0.5",         # Start distance
            "3.0",         # End distance
            "0.01",        # Step size
            $outFile       # Output filename
        )

        # Put them in the script and run
        $Inputs | python hf.py
        Write-Host "`nFinished writing to $outFile`n" -ForegroundColor Green
    }
}
Write-Host "All batch calculations complete." -ForegroundColor Cyan