#!/bin/bash -l

#SBATCH --qos=normal
#SBATCH --partition=shared-hopper-devkit
#SBATCH --job-name=srm_sens
#SBATCH --time=10:00:00

# Array of dt values
dt_values=(0.00005)

# Loop through each dt value
for dt in "${dt_values[@]}"
do
    dt_str="${dt/./p}"
    echo "Running with dt = $dt"
    srun -n 1 -G 1 julia src/single_resonance/magnus_tau_dep.jl $dt > "output_dt_${dt_str}.txt"
done