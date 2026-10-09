module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions, AtomicAndPhysicalConstants, LinearAlgebra
using Beamlines: isactive, deval, unsafe_getparams, isnullspecies
using BeamTracking: R_to_E, R_to_beta_gamma, R_to_gamma, R_to_pc, R_to_v,
                    beta_gamma_to_v, E_to_R, E_to_v,
                    @makekernel, Coords, make_kernel_call, KernelCall, KernelChain, push, TimeDependentParam, RefState, 
                    launch!, AbstractSymplectic, rot_quaternion, inv_rot_quaternion, atan2, 
                    get_N_particle, mean_and_cov, ibs_integrals, remake, push_transforms_in, push_transforms_out
                    
import BeamTracking: track!

function track!(
  bunch::Bunch, 
  ele::LineElement;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  _p_over_q_ref=nothing,
  context=(haskey(getfield(ele, :pdict), BeamlineParams) ? ele.beamline.context : Beamlines.NULL_CONTEXT),
  batch_start=1,
  kwargs...
)
  if isnothing(_p_over_q_ref)
    p_over_q_ref = ele.p_over_q_ref
  else
    p_over_q_ref = _p_over_q_ref
  end
  coords = bunch.coords
  @noinline _track!(coords, bunch, ele, context, p_over_q_ref, ele.tracking_method, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start; kwargs...)
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
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end
  __, p_over_q_ref = check_bl_bunch!(bunch, bl)
  context = bl.context
  
  for ele in bl.line
    track!(bunch, ele; context, scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start, kwargs...)
  end

  return bunch
end

function track!(
  bunch::Bunch, 
  branch::Branch; 
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  batch_start=1,
  kwargs...
)
  # Empty Beamlines are skipped: they contain no elements to track through
  bls = filter(bl -> length(bl.line) != 0, branch.beamlines)
  if length(bls) == 0
    return bunch
  end
  # Set the bunch species and reference energy from the start of the Branch if not set
  check_bl_bunch!(bunch, first(bls))

  for (i, bl) in enumerate(bls)
    if i > 1
      step_reference_energy!(bunch, bls[i-1], bl; kwargs...)
    end
    track!(bunch, bl; scalar_params, ramp_particle_energy_without_rf, ramp_update_each_particle, rf_on, batch_start, kwargs...)
  end

  return bunch
end

"""
    step_reference_energy!(bunch::Bunch, bl_up::Beamline, bl_down::Beamline; kwargs...)

Shifts the reference energy of `bunch` by the step in reference energy from Beamline `bl_up` 
to the directly-downstream Beamline `bl_down`, with both reference energies evaluated at 
`bunch.t_ref`. The phase space coordinates are rescaled to the new reference energy keeping 
the particle energies unchanged, independent of `ramp_particle_energy_without_rf`, since a step 
between Beamlines is a change in reference only. Any remaining difference between `bunch.p_over_q_ref` 
and the reference energy of `bl_down` (from ramping) is handled by element tracking. 
`kwargs` are passed to `launch!`.
"""
function step_reference_energy!(bunch::Bunch, bl_up::Beamline, bl_down::Beamline; kwargs...)
  p_over_q_up = Beamlines.trygetproperty(bl_up, :p_over_q_ref)
  p_over_q_down = Beamlines.trygetproperty(bl_down, :p_over_q_ref)
  if p_over_q_up isa Beamlines.GetError || p_over_q_down isa Beamlines.GetError
    return bunch
  end

  t_ref = bunch.t_ref
  p_over_q_up = p_over_q_up isa TimeDependentParam ? p_over_q_up(t_ref) : p_over_q_up
  p_over_q_down = p_over_q_down isa TimeDependentParam ? p_over_q_down(t_ref) : p_over_q_down
  dp_over_q = p_over_q_down - p_over_q_up
  if dp_over_q == 0
    return bunch
  end

  p_over_q_old = bunch.p_over_q_ref
  p_over_q_new = p_over_q_old + dp_over_q
  beta_gamma_old = R_to_beta_gamma(bunch.species, p_over_q_old)
  dbeta_gamma = R_to_beta_gamma(bunch.species, p_over_q_new) - beta_gamma_old
  T = eltype(bunch.coords.v)
  kcall = make_kernel_call(BeamTracking.reference_momentum_shift!, 
            (BeamTracking.num_lower(T, beta_gamma_old), BeamTracking.num_lower(T, dbeta_gamma), Val{true}()))
  @noinline launch!(bunch.coords, kcall; kwargs...)
  bunch.p_over_q_ref = p_over_q_new
  return bunch
end

include("utils_bl.jl")
include("unpack_bl.jl")
include("scibmadstandard_bl.jl")
include("exact_bl.jl")
include("symplectic_bl.jl")
include("sagan_cavity_bl.jl")
include("general_bl.jl")
include("rungekutta.jl")

end