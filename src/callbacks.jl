#=
Callbacks must have the signature

i, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in

where transforms_out and transforms_in are KernelChain of kernels with signatures
i, coords, cur_s

that must be applied in order to transform coords from the body frame to global frame (out)
and back from the global frame to the body frame (in)

Users can then decide in their callbacks if they 

=#

# For no callbacks, it doesn't matter what L, last_ds_step, and last_g 
# are since the callbacks will never be called.
process_callbacks(callbacks::Tuple{}, chain::T, callback_chain::S) where {T,S} = callbacks, 0, 0, 0 
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
@generated function process_callbacks(callbacks::C, chain::T, callback_chain::S) where {C, T, S}
  # Cases: Misalignments, Implicit, Ramping all require transformation back to global standard coordinates
  N = length(T.parameters)
  
  # First check if ramp_update_each_particle, in which case the first KernelCall will 
  # be a time-dependent reference_momentum_shift! that we must undo. This one is 
  # hard bc we need to know t_initial to undo it. For now, just don't support it.
  if first(T.parameters) <: KernelCall{typeof(reference_momentum_shift!),Tuple{<:Any,TimeFunction}}
    error("Callbacks with ramp_update_each_particle = true are not supported yet.")
  end

  # Remove apertures:
  if first(T.parameters) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(process_callbacks(callbacks, Base.tail(chain), callback_chain))
  elseif last(T.parameters) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(process_callbacks(callbacks, Base.front(chain), callback_chain))
  end

  # Misalignments:
  if first(T.parameters) <: KernelCall{typeof(track_coord_transform!)}

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
execute_callbacks(i, coords, cur_s, last_ds_step, last_g) = _execute_callbacks(i, coords.callbacks, coords, cur_s, last_ds_step, last_g)

@unroll function _execute_callbacks(i, callbacks, coords, cur_s, last_ds_step, last_g)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, last_ds_step, last_g)
  end
  return nothing
end
