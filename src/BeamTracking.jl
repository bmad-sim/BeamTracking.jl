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
      AtomicAndPhysicalConstants

using KernelAbstractions
using SIMD: SIMD
      
import GTPSA: sincu, sinhcu, normTPS
import Base: setproperty!

export Bunch, State, ParticleView, Time, TimeDependentParam
export LinearTracking, Linear, ExactTracking, Exact
export IntegrationTracking, SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick
export track!
export sincu, sinhcu, sincuc, expq, atan2, one_cos, one_cos_norm
export rot_quaternion, inv_rot_quaternion, quat_mul, quat_rotate
export Species

include("utils.jl")
include("types.jl")
include("time.jl")
include("kernel.jl")

include("helpers/track_transforms.jl")
include("helpers/coord_transforms.jl")
include("kernels/alignment_kernel.jl")
include("kernels/aperture_kernel.jl")
include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear
include("modules/IntegrationTracking.jl") #; TRACKING_METHOD(::LinearTracking) = SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick

# Empty tracking method to be imported+implemented by package extensions
function track! end

# --------------------------------------------------

# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")

end
