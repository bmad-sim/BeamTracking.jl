# BeamTracking

[![Build Status](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides universally polymorphic and portable, parallelizable routines for simulating charged particle beams both on the CPU and, using [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl), various GPU backends including NVIDIA CUDA, Apple Metal, Intel oneAPI, and AMD ROCm.

`BeamTracking.jl` is one package in the [SciBmad](https://github.com/bmad-sim/SciBmad.jl) ecosystem of accelerator physics software. 

All implemented integrators have been tested against the [Polymorphic Tracking Code (PTC)](https://cds.cern.ch/record/573082/files/sl-2002-044.pdf).

See [the SciBmad documentation](https://bmad-sim.github.io/SciBmad.jl/) for documentation.