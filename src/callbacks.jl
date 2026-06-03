
# For no callbacks, it doesn't matter what L, last_ds_step, and last_g 
# are since the callbacks will never be called.
process_callbacks(callbacks::Tuple{}, chain::T) where {T} = callbacks, 0, 0, 0 
#=
@generated function beval(f::T, t) where {T<:Tuple}
  N = length(T.parameters)
  # Use getfield with literal integer arguments
  exprs = [:(beval(Base.getfield(f, $i), t)) for i in 1:N]
  return :(tuple($(exprs...)))
end
=#
# This function looks at a Tuple of KernelCalls and returns 
# new set of callbacks that acts on coords including automatic 
# transformation back to 
@generated function process_callbacks(callbacks::C, chain::T) where {C, T}
  # Cases: Misalignments, Implicit, Ramping all require transformation back to global standard coordinates
  N = length(T.parameters)
  
  # First check if ramp_update_each_particle, in which case the first KernelCall will 
  # be a time-dependent reference_momentum_shift! that we must undo. This one is 
  # hard bc we need to know t_initial to undo it. For now, just don't support it.
  if first(T) <: KernelCall{typeof(reference_momentum_shift!),Tuple{<:Any,TimeFunction}}
    error("Callbacks with ramp_update_each_particle = true are not supported yet.")
  end

  # Remove apertures:
  if first(T) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(transform_callbacks(callbacks, Base.tail(chain)))
  elseif last(T) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(transform_callbacks(callbacks, Base.front(chain)))
  end

  # Misalignments:
  if first(T) <: KernelCall{typeof(track_coord_transform!)}

  end
  

  if N > 0 && first(T.parameters) <: KernelCall{typeof(reference_momentum_shift!),Tuple{<:Any,TimeFunction}}
    # Static check that everything is ok
    if last(T.parameters) <: KernelCall{typeof(reference_momentum_shift!),Tuple{TimeFunction,TimeFunction}}
      return :(__generic_kernel_ramp!(i, coords, chain, ref))
    else
      error("
        Kernels with time-dependent reference energies must start and end with `reference_momentum_shift!`,
        where the entering one transforms all particles to have individual reference momenta, and the exiting 
        transforms all particles to a uniform reference momentum at the end of the element equal to 
        p_over_q_ref(bunch.t_ref at end).
      ")
    end
  else
    return :(__generic_kernel_noramp!(i, coords, chain, ref))
  end
end

# Function to execute callbacks
execute_callbacks(coords, cur_s, last_ds_step, last_g) = _execute_callbacks(coords.callbacks, coords, cur_s, last_ds_step, last_g)

@unroll function _execute_callbacks(callbacks, coords, cur_s, last_ds_step, last_g)
  @unroll for callback in callbacks
    callback(coords, cur_s, last_ds_step, last_g)
  end
  return nothing
end
