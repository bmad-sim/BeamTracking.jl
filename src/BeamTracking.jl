module BeamTracking
using SciMLBase,
      RecursiveArrayTools,
      NaNMath,
      GTPSA,
      ReferenceFrameRotations,
      StaticArrays,
      StructArrays

include("aapc.jl")

export Bunch,       
       Coord, 
       Quaternion,
       Particle,
       
       MatrixKick,
       Linear,
       DiffEq,

       sr_gamma, 
       sr_gamma_m1,
       sr_beta,
       sr_pc,
       sr_ekin,
       sr_etot,
       brho,

       chargeof,
       massof,
       Species,

       sincu,
       sinhcu,

       track!





include("new_types.jl")

# Empty tracking method to be imported by submodules 
track!(bunch::Bunch, ::Nothing; work=nothing) = bunch

# --------------------------------------------------

include("utils.jl")
include("work.jl")

# Modules separated:
include("DiffEq/DiffEq.jl")
include("MatrixKick/MatrixKick.jl") 
include("Linear/Linear.jl")   
include("Misc/Misc.jl")



end
