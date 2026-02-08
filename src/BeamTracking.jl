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
      EnumX

using KernelAbstractions

import GTPSA: sincu, sinhcu

export Bunch, State, ParticleView, Time, TimeDependentParam, BatchParam
export Yoshida, Yoshida, MatrixKick, BendKick, SolenoidKick, DriftKick, Exact
export Fringe, SaganCavity, track!


include("utils/coord_transforms.jl")
include("utils/energy.jl")
include("utils/math_simd.jl")
include("utils/quaternions.jl")
include("utils/z_to_time.jl")

include("types.jl")
include("time.jl")
include("batch.jl")
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
include("kernels/ramp_P0.jl")
include("kernels/rfcavity_kick.jl")
include("kernels/sagan_cavity.jl")
include("kernels/solenoid_kick.jl")
include("kernels/spin.jl")
include("kernels/transforms.jl")
include("kernels/yoshida.jl")

# Empty tracking method to be imported+implemented by package extensions
function track! end

end
