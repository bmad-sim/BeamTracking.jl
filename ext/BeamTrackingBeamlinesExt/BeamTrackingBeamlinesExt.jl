module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions, AtomicAndPhysicalConstants
using Beamlines: isactive, deval, unsafe_getparams, isnullspecies
using BeamTracking: R_to_E, R_to_beta_gamma, R_to_gamma, R_to_pc, R_to_v, 
                    beta_gamma_to_v, E_to_R, E_to_v,
                    @makekernel, Coords, KernelCall, KernelChain, push, TimeDependentParam, RefState, 
                    launch!, AbstractYoshida, rot_quaternion, inv_rot_quaternion, atan2, 
                    get_N_particle
                    
import BeamTracking: track!

include("utils.jl")

function track!(
  bunch::Bunch, 
  ele::LineElement;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  kwargs...
)
  coords = bunch.coords
  @noinline _track!(coords, bunch, ele, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf; kwargs...)
  return bunch
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end
  check_bl_bunch!(bl, bunch)
  
  for ele in bl.line
    track!(bunch, ele; scalar_params, ramp_particle_energy_without_rf, kwargs...)
  end

  return bunch
end

include("rfcavity_bl.jl")
include("unpack.jl")
include("scibmadstandard.jl")
include("exact.jl")
include("yoshida.jl")
include("sagan_cavity_bl.jl")
include("general_bl.jl")

end