module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      VectorizationBase,
      EnumX,
      Unrolled,
      MacroTools,
      Adapt,
      Accessors

using KernelAbstractions
      
import GTPSA: sincu, sinhcu
import Base: setproperty!

export Bunch, Species, State, ParticleView, ELECTRON, POSITRON, PROTON, ANTIPROTON, sincu, sinhcu, sincuc
export LinearTracking, Linear
export ExactTracking, Exact
export track!, track_chain

include("utils.jl")
include("types.jl")
include("kernel.jl")



include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear


# Empty tracking method to be imported+implemented by package extensions
function track! end

# To get a track_chain
function track_chain end

function MAX_TEMPS end
# --------------------------------------------------


# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")


end
