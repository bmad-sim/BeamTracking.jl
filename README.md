# BeamTracking

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/BeamTracking.jl/)
[![Build Status](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides universally polymorphic and portable, parallelizable routines for simulating charged particle beams both on the CPU and, using [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl), various GPU backends including NVIDIA CUDA, Apple Metal, Intel oneAPI, and AMD ROCm.

`BeamTracking.jl` is one package in the [SciBmad](https://github.com/bmad-sim/SciBmad.jl) ecosystem of accelerator physics software. 

Currently, the following integrators have been implemented and tested against the [Polymorphic Tracking Code (PTC)](https://cds.cern.ch/record/573082/files/sl-2002-044.pdf):

- Yoshida's second, fourth, sixth, and eighth-order symplectic methods for the following splits:
    - Drift-kick
    - Thick-kick (including quadrupole-kick, solenoid-kick, bend-kick)
- Exact solenoid, exact bend, patch


For examples of how to use the package, we currently refer users to the [`SciBmad.jl`](https://github.com/bmad-sim/SciBmad.jl) package while the documentation is still in development.