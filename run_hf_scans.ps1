# run_hf_scans.ps1

# Define our molecules and the input key your script expects
$molecules = @(
    #@{ Key="a"; Name="H2" },
    @{ Key="b"; Name="HeH+" }
)

# Loop through both molecules
foreach ($mol in $molecules) {
    
    # Loop through basis sets n = 2, 3, 4, 5, 6
    for ($n = 2; $n -le 6; $n++) {
        
        # Create a unique output file name for this run
        $outFile = "$($mol.Name)_STO-$($n)G.csv"
        
        Write-Host "======================================================" -ForegroundColor Cyan
        Write-Host " Running $($mol.Name) with STO-$($n)G basis set..." -ForegroundColor Cyan
        Write-Host "======================================================" -ForegroundColor Cyan

        # Create an array of strings in the exact order your script asks for them
        $simulatedInputs = @(
            $mol.Key,      # Select Molecule (a or b)
            $n.ToString(), # Basis set (2-6)
            "0.5",         # Start distance
            "3.0",         # End distance
            "0.01",        # Step size
            $outFile       # Output filename
        )

        # Pipe the array into your python script
        $simulatedInputs | python hf.py
        
        Write-Host "`nFinished writing to $outFile`n" -ForegroundColor Green
    }
}

Write-Host "All batch calculations complete." -ForegroundColor Cyan