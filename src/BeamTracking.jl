"""
    BeamTracking

A high-performance particle beam tracking package for accelerator physics simulations.
Currently provides both linear, exact, field tracking, and Runge-Kutta tracking methods.
"""
module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      VectorizationBase

using KernelAbstractions
      
import GTPSA: sincu, sinhcu
import Base: setproperty!

export Bunch, Species, ParticleView, ELECTRON, POSITRON, PROTON, ANTIPROTON, sincu, sinhcu, sincuc
export LinearTracking, Linear
export ExactTracking, Exact
export FieldTracking, Field
export RungeKuttaTracking, RungeKutta
export track!

include("utils.jl")
include("kernel.jl")
include("types.jl")


include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear
include("modules/FieldTracking.jl") #; TRACKING_METHOD(::FieldTracking) = Field
include("modules/RungeKuttaTracking.jl") #; TRACKING_METHOD(::RungeKuttaTracking) = RungeKutta

# Empty tracking method to be imported+implemented by package extensions
function track! end

function MAX_TEMPS end
# --------------------------------------------------


# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")


end
