# BeamTracking

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/BeamTracking.jl/)
[![Build Status](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)

A high-performance charged particle beam tracking package that provides tracking methods for accelerator physics simulations. Features include:

- Universal polymorphism and portability across CPU and GPU backends (NVIDIA CUDA, Apple Metal, Intel oneAPI, AMD ROCm)
- Parallel processing with SIMD, multi-threading, and [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl)
- Extensive collection of tracking methods

To develop this package:

```julia
import Pkg;
Pkg.develop(url="https://github.com/bmad-sim/BeamTracking.jl.git"); # This package! Replace bmad-sim with your username if working on a fork
```

If working on your own fork, replace `bmad-sim` in the above `develop` url with your Github username.

In your `~/.julia/dev/` directory, you will now see the directory `BeamTracking`. This is the Github repo where you can do your work and push changes.

See the [development documentation](https://bmad-sim.github.io/BeamTracking.jl/dev/) for more details.

## Tracking Methods

1. **Linear Tracking** (`src/modules/LinearTracking.jl`)
2. **Exact Tracking** (`src/modules/ExactTracking.jl`)
