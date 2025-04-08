module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      VectorizationBase,
      Polyester

import GTPSA: sincu, sinhcu
import Base: setproperty!

export Bunch, Species, ParticleView, ELECTRON, POSITRON, PROTON, ANTIPROTON, sincu, sinhcu
export LinearTracking, Linear
export ExactTracking, Exact
export track!

include("utils.jl")
include("kernel.jl")
include("types.jl")

include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear

# Empty tracking method to be imported+implemented by package extensions
function track! end

function MAX_TEMPS end
# --------------------------------------------------


# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")


end
