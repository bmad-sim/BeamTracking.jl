#!/bin/bash

# Loop over process counts (-n) and thread counts (-t)
for n in 2 4 8; do
    for t in 32 64 128 256; do
        echo "Running with -n $n and -t $t"
        # srun -n $n --cpu-bind=cores julia -t $t mpi_example_no_plot.jl 100000
        mpiexecjl -n $n julia -t $t examples/mpi/mpi_example.jl 100000
    done
done