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
  VectorizationBase,
  Unrolled,
  MacroTools,
  Adapt,
  Accessors,
  SIMDMathFunctions

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
export FieldTracking, Field
export RungeKuttaTracking, RungeKutta
export IntegrationTracking, SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick
export track!

include("utils.jl")
include("types.jl")
include("kernel.jl")



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
