#!/bin/bash -l

# Array of dt values
nu_values=(62.47 62.48 62.49)

# Loop through each dt value
for nu in "${nu_values[@]}"
do
    nu_str="${nu}"
    echo "Running with nu = $nu"
    julia src/single_resonance/magnus_tau_dep.jl $nu > "output_nu_${nu_str}.txt"
done