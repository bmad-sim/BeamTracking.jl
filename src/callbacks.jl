#=
Callbacks must have the signature

i, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in

where transforms_out and transforms_in are KernelChains of kernels with signatures
i, coords, cur_s

that must be applied in order to transform coords from the body frame to global frame (out)
and back from the global frame to the body frame (in)

Users can then decide in their callbacks if they 
Now for each callback in the callback chain, we need to construct closures
where transforms_in and transforms_out are provided. e.g.
mycallback(i, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in)
needs to be closed:
let transforms_out=..., tranforms_in=...
  (i, coords, cur_s, last_ds_step, last_g) -> mycallback(i, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in)
end
This means that we should do this in a two step process: 
1. construct transforms_out and transforms_in KernelChains
2. Construct a single callback that evaluates the transforms_out and transforms_in KernelChains at 
the given cur_s, and then evaluates all callbacks using those 
=#
@generated function construct_transforms(chain::S, transforms_out::T, transforms_in::U) where {S,T,U}
  N = length(S.parameters)
  if N == 0
    return :(transforms_out, transforms_in)
  end

  # First check if ramp_update_each_particle, in which case the first KernelCall will 
  # be a time-dependent reference_momentum_shift! that we must undo. This one is 
  # hard bc we need to know t_initial to undo it. For now, just don't support it.
  if first(S.parameters) <: KernelCall{typeof(reference_momentum_shift!),Tuple{<:Any,TimeFunction}}
    error("Callbacks with ramp_update_each_particle = true are not supported yet.")
  end

  # Remove apertures:
  if first(S.parameters) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(construct_transforms(Base.tail(chain), transforms_out, transforms_in))
  elseif last(S.parameters) <: KernelCall{<:Union{typeof(track_aperture_elliptical!),typeof(track_aperture_rectangular!)}}
    return :(construct_transforms(Base.front(chain), transforms_out, transforms_in))
  end

  # Misalignments:
  # Straight:
  if first(S.parameters) <: KernelCall{typeof(track_alignment_straight_at_s!)} && 
    last(S.parameters) <: KernelCall{typeof(track_alignment_straight_at_s!)}
    return quote
      #=
      @assert x_off == chain[end].args[1]
      @assert y_off == chain[ends].args[2]
      @assert z_off == chain[end].args[3]
      @assert x_rot == chain[end].args[4]
      @assert y_rot == chain[end].args[5]
      @assert tilt  == chain[end].args[6]
      @assert L     == chain[end].args[8]
      =#
      let x_off = chain[1].args[1], y_off = chain[1].args[2], z_off = chain[1].args[3], x_rot = chain[1].args[4], y_rot = chain[1].args[5], tilt  = chain[1].args[6], ele_orient = chain[1].args[7], L = chain[1].args[8]
        tin  = (i, coords, cur_s) -> track_alignment_straight_at_s!(i, coords, x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, cur_s, Val{true}())
        tout = (i, coords, cur_s) -> track_alignment_straight_at_s!(i, coords, x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, cur_s, Val{false}())
      end
      # Note that this KernelChain can't just be launch!-ed yet, because cur_s is not provided
      # execute_callbacks will handle this by modifying this KernelChain so cur_s is included in 
      # args before the callback is executed, so that in the callback the user can launch! it.
      transforms_out = push(transforms_out, KernelCall(tout, ()))
      transforms_in = push(transforms_in, KernelCall(tin, ()))
      construct_transforms(Base.front(Base.tail(chain)), transforms_out, transforms_in)
    end
  end 
  
  return :(transforms_out, transforms_in)
end

function construct_main_callback(callbacks, _transforms_out, _transforms_in)
  let _transforms_in=_transforms_in, _transforms_out=_transforms_out, callbacks=callbacks
    return (i, coords, cur_s, last_ds_step, last_g) -> begin
      transforms_out = KernelChain(_add_s_to_chain(cur_s, _transforms_out.chain), _transforms_out.ref)
      transforms_in = KernelChain(_add_s_to_chain(cur_s, _transforms_in.chain), _transforms_in.ref)
      _execute_callbacks(i, callbacks, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in)
      return nothing
    end
  end
end

@unroll function _add_s_to_chain(s, chain)
  i = 0
  @unroll for kcalli in chain
    i += 1
    if kcalli.kernel != blank_kernel!
      @reset chain[i] = KernelCall(kcalli, (s))
    end
  end
  return chain
end

@unroll function _execute_callbacks(i, callbacks, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, last_ds_step, last_g, transforms_out, transforms_in)
  end
  return nothing
end

# Function to execute main_callbac
# after addition of transform_out and transform_in, this reall just executes the main callback
execute_callbacks(i, coords, cur_s, last_ds_step, last_g) = __execute_callbacks(i, coords.callbacks, coords, cur_s, last_ds_step, last_g)

@unroll function __execute_callbacks(i, callbacks, coords, cur_s, last_ds_step, last_g)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, last_ds_step, last_g)
  end
  return nothing
end

#=

THIS IS NOT A NICE SOLUTION! We should explore more general solutions rather than 
KernelChain introspection. Nonetheless, it works for now.

=#
@generated function find_L_ds_step_and_g(chain::C) where {C}
  if length(C.parameters) == 0
    return :(0, 0, 0)
  end

  yoshida_kcall_idx = findfirst(x->x <: typeof(order_two_integrator!) || x <: typeof(order_four_integrator!) ||
        x <: typeof(order_six_integrator!) || x <: typeof(order_eight_integrator!), C.parameters)
  if !isnothing(yoshida_kcall_idx)
    return quote
      ker = chain[$yoshida_kcall_idx].args[1]
      params = chain[$yoshida_kcall_idx].args[2]
      L = chain[$yoshida_kcall_idx].args[9]
      ds_step = chain[$yoshida_kcall_idx].args[4]
      g = compute_g(ker, params)
      return L, ds_step, g
    end
  end

  thin_kcall_idx = findfirst(x->x <: typeof(integrate_thin!), C.params)
  if !isnothing(thin_kcall_idx)
    return :(0, 0, 0)
  end


  if any(x -> x <: typeof(order_two_integrator!) || x <: typeof(order_four_integrator!) ||
        x <: typeof(order_six_integrator!) || x <: typeof(order_eight_integrator!))
    kcall = findfirst(x->x)
    # ds_step is 
  end
end