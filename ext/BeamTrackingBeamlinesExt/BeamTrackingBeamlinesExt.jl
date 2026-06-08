module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions, AtomicAndPhysicalConstants, LinearAlgebra
using Beamlines: isactive, deval, unsafe_getparams, isnullspecies
using BeamTracking: R_to_E, R_to_beta_gamma, R_to_gamma, R_to_pc, R_to_v, 
                    beta_gamma_to_v, E_to_R, E_to_v,
                    @makekernel, Coords, make_kernel_call, KernelCall, KernelChain, push, TimeDependentParam, RefState, 
                    launch!, AbstractYoshida, rot_quaternion, inv_rot_quaternion, atan2, 
                    get_N_particle, mean_and_cov, ibs_integrals, remake
                    
import BeamTracking: track!

function track!(
  bunch::Bunch, 
  ele::LineElement;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  _p_over_q_ref=nothing,
  kwargs...
)
  if isnothing(_p_over_q_ref)
    p_over_q_ref = ele.p_over_q_ref
  else
    p_over_q_ref = _p_over_q_ref
  end
  coords = bunch.coords
  @noinline _track!(coords, bunch, ele,  p_over_q_ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle; kwargs...)
  return bunch
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end
  __, p_over_q_ref = check_bl_bunch!(bunch, bl)
  
  for ele in bl.line
    track!(bunch, ele; scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, kwargs...)
  end

  return bunch
end

include("utils_bl.jl")
include("unpack_bl.jl")
include("scibmadstandard_bl.jl")
include("exact_bl.jl")
include("yoshida_bl.jl")
include("sagan_cavity_bl.jl")
include("general_bl.jl")

end