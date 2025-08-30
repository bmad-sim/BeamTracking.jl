module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      SIMDMathFunctions,
      VectorizationBase,
      Unrolled,
      MacroTools,
      Adapt,
      Accessors

using KernelAbstractions
      
import GTPSA: sincu, sinhcu, normTPS
import Base: setproperty!
import SIMD: Vec

# Put AtomicAndPhysicalConstants in a box for now for safety
include("Constants.jl")
using .Constants: Constants, Species, massof, chargeof, nameof, C_LIGHT, isnullspecies
export Species

export Bunch, State, ParticleView, sincu, sinhcu, sincuc, expq, quat_mul, atan2
export LinearTracking, Linear
export ExactTracking, Exact
export IntegrationTracking, SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick
export track!
export rot_quaternion, inv_rot_quaternion

include("utils.jl")
include("types.jl")
include("kernel.jl")
include("kernels/alignment.jl")
include("kernels/aperture.jl")
include("kernels/coord_rotation.jl")

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
