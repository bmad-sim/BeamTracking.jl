# BeamTracking

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/BeamTracking.jl/stable/)
[![Build Status](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)

A high-performance particle beam tracking package written in Julia, designed for accelerator physics simulations.

## Overview

BeamTracking is a specialized package for tracking particle beams through accelerator elements. It provides various tracking methods and memory layouts optimized for performance.

## Installation

To install the package:

```julia
using Pkg
Pkg.add("BeamTracking")
```

To develop this package:

```julia
import Pkg;
Pkg.develop(url="https://github.com/bmad-sim/BeamTracking.jl.git"); # This package! Replace bmad-sim with your username if working on a fork
```

In your `~/.julia/dev/` directory, you will now see the directory `BeamTracking`. This is the Github repo where you can do your work and push changes.

See the [development documentation](https://bmad-sim.github.io/BeamTracking.jl/dev/) for more details.

## Core Components

### Data Structures

- `Bunch`: The main data structure representing a particle bunch
  - Supports both AoS and SoA memory layouts
  - Contains species information and reference magnetic rigidity (Brho_ref)
  - Stores particle coordinates in a matrix format

### Tracking Methods

1. **Linear Tracking** (`LinearTracking.jl`)
2. **Exact Tracking** (`ExactTracking.jl`)

### Performance Features

- SIMD vectorization support
- Automatic multithreading for large particle numbers
- Optimized memory layouts (AoS/SoA)
- Efficient kernel launching system
- GPU acceleration support. Compatible with KernelAbstractions.jl
- Seamless CPU/GPU code sharing through unified kernel interface
