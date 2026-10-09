module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      SIMDMathFunctions,
      HostCPUFeatures,
      Unrolled,
      MacroTools,
      Adapt,
      Accessors,
      SpecialFunctions,
      AtomicAndPhysicalConstants,
      Random,
      EnumX,
      Statistics,
      LinearAlgebra,
      TPSAInterface,
      ForwardDiff

using KernelAbstractions

import GTPSA: sincu, sinhcu

export Species
export Bunch, State, ParticleView, Time, TimeDependentParam, BatchParam
export Symplectic, MatrixKick, BendKick, SolenoidKick, DriftKick, Exact, RungeKutta
export Fringe, SaganCavity, track!

include("utils/float32.jl")
include("utils/coord_transforms.jl")
include("utils/energy.jl")
include("utils/math_simd.jl")
include("utils/quaternions.jl")
include("utils/z_to_time.jl")
include("utils/beam_statistics.jl")
include("utils/ibs_integrals.jl")
include("utils/random.jl")

include("types.jl")
include("time.jl")
include("batch.jl")
include("callbacks.jl")
include("kernels/ramp_P0.jl")
include("kernel.jl")
include("tracking_methods.jl")

include("kernels/kernel_utils.jl")
include("kernels/alignment.jl")
include("kernels/aperture.jl")
include("kernels/bend_kick.jl")
include("kernels/drift_kick.jl")
include("kernels/map.jl")
include("kernels/multipole.jl")
include("kernels/patch.jl")
include("kernels/quadrupole_kick.jl")
include("kernels/radiation.jl")
include("kernels/rfcavity_kick.jl")
include("kernels/sagan_cavity.jl")
include("kernels/solenoid_kick.jl")
include("kernels/spin.jl")
include("kernels/transforms.jl")
include("kernels/integrators.jl")
include("kernels/ibs_kick.jl")
include("kernels/implicit.jl")
include("kernels/thin.jl")
include("kernels/fringe.jl")
include("kernels/elsep.jl")

include("fields.jl")
include("kernels/runge_kutta.jl")

include("utils/find_stuff.jl")

# 
"""
    track!(bunch::Bunch, ele::LineElement; kwargs...) -> bunch
    track!(bunch::Bunch, bl::Beamline; kwargs...) -> bunch
    track!(bunch::Bunch, branch::Branch; kwargs...) -> bunch

Tracks the particles in `bunch` in place through a single `LineElement` or all the
elements of a `Beamline`, or through all the `Beamline`s of a `Branch`.

## Reference species and energy

When tracking through a `Beamline` or `Branch`: If `bunch.species` is not set, it is set to
the reference species of the `Beamline` (for a `Branch`, the first non-empty `Beamline`). If it
is set, it must equal the reference species of every `Beamline` tracked through. Similarly, if
`bunch.p_over_q_ref` is `NaN`, it is set to the reference `p_over_q_ref` of the `Beamline` (or
the first non-empty `Beamline` of the `Branch`). When tracking through a single `LineElement`,
these are not set from the element, so `bunch.species` and `bunch.p_over_q_ref` should already
be set.

## Keyword arguments

- `scalar_params::Bool = false`: If `true`, element parameters are converted to regular
  number types (see `Beamlines.scalarize`) before tracking. Useful when parameters carry
  e.g. `GTPSA` or `ForwardDiff` types and their derivatives are not needed.
- `ramp_particle_energy_without_rf::Bool = false`: When the reference energy of an element
  differs from `bunch.p_over_q_ref` (e.g. when the reference energy is ramped with a
  `TimeDependentParam`), if `true` the particle energies are scaled with the reference energy
  so the energy deviations are unchanged. If `false`, the particle energies are unchanged and
  only the coordinates are rescaled.
- `ramp_update_each_particle::Bool = false`: If `true` and the reference energy is a
  `TimeDependentParam`, the reference energy is evaluated at the time of each particle
  separately. If `false`, it is evaluated at `bunch.t_ref` at the element entrance.
- `rf_on::Bool = true`: If `false`, `RFParams` are ignored so RF cavities are tracked as
  drifts (with any other parameters of the element).
- `batch_start = 1`: Index offset into `BatchParam`s: particle `i` uses the
  `mod1(i + batch_start - 1, N_batch)`-th value of each `BatchParam`.
- `context`: (`LineElement` only) The `Context` holding the control variables used to
  evaluate deferred expressions. Defaults to that of the element's `Beamline`, if any. For a
  `Beamline` or `Branch` the `Context` of the `Beamline` is used.

Remaining keyword arguments are passed to the kernel launcher. These are:
- `use_KA::Bool`: Use `KernelAbstractions` to launch the kernels. Defaults to `true` unless
  the coordinates are on the CPU and `groupsize` is not set.
- `use_explicit_SIMD::Bool = !use_KA`: Use explicit SIMD vectorization on the CPU.
- `use_cpu_multithreading::Bool = false`: Multithread over the particles on the CPU.
- `groupsize::Union{Nothing,Integer} = nothing`: `KernelAbstractions` workgroup size.

## Example
```julia
using Beamlines, BeamTracking
qf = Quadrupole(Kn1=0.36, L=0.5)
d = Drift(L=1)
qd = Quadrupole(Kn1=-0.36, L=0.5)
fodo = Beamline([qf, d, qd, d], species_ref=Species("electron"), E_ref=18e9)

bunch = Bunch(zeros(10, 6))
bunch.v[:, BeamTracking.XI] .= 1e-3
track!(bunch, fodo)  # bunch.species and bunch.p_over_q_ref are set from fodo
```
"""
function track! end

end
