module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions, AtomicAndPhysicalConstants, LinearAlgebra
using Beamlines: isactive, deval, R_to_E, R_to_beta_gamma, R_to_pc, R_to_v, E_to_R, E_to_v, LatElement, AbstractSymplectic
using BeamTracking: beta_gamma_to_v, make_kernel_call, KernelChain, push, TimeDependentParam, RefState, 
                    launch!, rot_quaternion, inv_rot_quaternion, atan2, exact_wedge!,
                    mean_and_cov, ibs_integrals, push_transforms_in, push_transforms_out, @set
                    
import BeamTracking: track!

function track!( 
  bunch::Bunch, 
  ele::LatElement;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  context=(haskey(getfield(ele, :pdict), BeamlineParams) ? ele.beamline.context : Beamlines.NULL_CONTEXT),
  batch_start=1,
  notify::Bool=true,
  kwargs...
)
  if haskey(getfield(ele, :pdict), BeamlineParams)
    bl = ele.beamline
    ref = resolve_reference!(bunch, Beamlines.trygetproperty(bl, :species_ref), Beamlines.trygetproperty(bl, :p_over_q_ref), notify)
  else
    # Standalone element: everything from the Bunch
    ref = resolve_reference!(bunch, nothing, nothing, notify)
  end
  @noinline _track!(bunch, ele, context, ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start; kwargs...)
  return bunch
end


function track!(
  bunch::Bunch, 
  bl::Beamline; 
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  batch_start=1,
  notify::Bool=true,
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end
  # Resolved once for the whole Beamline
  ref = resolve_reference!(bunch, Beamlines.trygetproperty(bl, :species_ref), Beamlines.trygetproperty(bl, :p_over_q_ref), notify)
  context = bl.context

  for ele in bl.line
    @noinline _track!(bunch, ele, context, ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start; kwargs...)
  end

  return bunch
end

include("utils_bl.jl")
include("unpack_bl.jl")
include("general_bl.jl")
include("symplectic_bl.jl")
include("driftkick_bl.jl")
include("matrixkick_bl.jl")
include("solenoidkick_bl.jl")
include("bendkick_bl.jl")
include("implicit_bl.jl")
include("exact_bl.jl")
include("rk4_bl.jl")
include("sagancavity_bl.jl")

end