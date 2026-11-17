#=
Example callback:

function mycallback(i, coords, p::TrackPoint)
  # Transform out of body frame:
  transforms_out!(i, coords, p)

  # Do stuff in global frame (with p.s, p.t_ref, p.step, ...)

  # Transform back into body frame:
  transforms_in!(i, coords, p)
end

=#

"""
    TrackPoint

The point of the tracking at which the callbacks are executed, passed to each callback as
`callback(i, coords, p::TrackPoint)`. Its properties are:

- `s`: s from the beginning of the `Beamline` (when tracking a single element, continued from the
  position of the `Bunch`: from the beginning of the element for a new `Bunch`)
- `t_ref`: reference time
- `step`: integration steps from the beginning of the `Beamline` (`n_slices*n_steps` per element,
  the last one at the exit of the element; when tracking a single element, continued from the
  position of the `Bunch` as `s`)
- `beta_gamma_ref`: reference energy (beta*gamma)
- `last_ds_step`: length of the last integration step
- `last_g`: curvature `(gx, gy)` of the last integration step

`transforms_out!(i, coords, p)` and `transforms_in!(i, coords, p)` transform particle `i` out of 
and back into the element body frame (e.g. of a misaligned element).

The transforms are only computed when used.
"""
struct TrackPoint{S,SE,T,D,B,L,G,TOUT,TIN}
  s::S                # from the Beamline entrance
  s_ele::SE           # from the element entrance (for the transforms)
  t_enter::T          # reference time at the launch entrance
  dt::D               # reference time since then
  step::Int
  beta_gamma_ref::B
  last_ds_step::L
  last_g::G
  batch_start::Int
  transforms_out::TOUT
  transforms_in::TIN
end

Base.getproperty(p::TrackPoint, key::Symbol) =
  key === :t_ref ? getfield(p, :t_enter) + getfield(p, :dt) : getfield(p, key)

Base.propertynames(::TrackPoint) = (:s, :t_ref, :step, :beta_gamma_ref, :last_ds_step, :last_g)

# Note: in keeping consistency with the rest of TimeDependentParam, transforms with time-dependent params 
# (e.g. misalignment with time-dependent x_offset) only evaluated at entrance of element. This is actually 
# needed to ensure everything is consistent, i.e., transform at exit consistent with rest of bunch. 
# perhaps we will want to refactor this if we want TimeDependentParam to evolve inside of elements, 
# but much lower priority.
transforms_out!(i, coords, p::TrackPoint) = _transform!(i, coords, p, getfield(p, :transforms_out))
transforms_in!(i, coords, p::TrackPoint) = _transform!(i, coords, p, getfield(p, :transforms_in))

function _transform!(i, coords, p, transforms)
  t_enter = getfield(p, :t_enter)
  args = _evaluate_transforms_args(i, coords, getfield(p, :batch_start), transforms, t_enter,
                                   getfield(p, :beta_gamma_ref))
  __merge_transforms(args, i, coords, getfield(p, :s_ele), p.t_ref)
  return nothing
end

# Constructs a single callback with signature (i, coords, cur_s, cur_dt_ref, cur_step) that executes 
# all callbacks in the callbacks tuple with the TrackPoint. cur_s and cur_step are from the launch 
# entrance: the callbacks get them from the Beamline entrance (s_enter, the s of the launch 
# entrance, and step0, the steps before it, are added), and the transforms get s from the element 
# entrance (s0, the s of the launch entrance within the element, nonzero for a slice of an element, 
# is added). cur_dt_ref is the reference time from the launch entrance.
function construct_main_callback(coords, batch_start, _transforms_out, _transforms_in, t_ref_transform, beta_gamma_ref_transform, ds_step, g, s0=0, s_enter=s0, step0=0)
  let batch_start = batch_start, _transforms_in = _transforms_in, _transforms_out = _transforms_out, t_ref_transform=t_ref_transform, beta_gamma_ref_transform=beta_gamma_ref_transform, callbacks=coords.callbacks, ds_step=ds_step, g=g, s0=s0, s_enter=s_enter, step0=step0
    return ((i, coords, cur_s, cur_dt_ref, cur_step) -> begin
      p = TrackPoint(s_enter+cur_s, s0+cur_s, t_ref_transform, cur_dt_ref, step0+cur_step, beta_gamma_ref_transform, ds_step, g, batch_start, _transforms_out, _transforms_in)
      _execute_callbacks_with_point(i, callbacks, coords, p)
      return nothing
    end,)
  end
end

@unroll function _evaluate_transforms_args(i, coords, batch_start, transforms, t_ref_transform, beta_gamma_ref_transform)
  idx = 1
  @unroll for transform in transforms
    bargs = process_batch_args(i, transform.args, batch_start)
    args = process_time_args(i, coords, bargs, t_ref_transform, beta_gamma_ref_transform)
    @reset transforms[idx].args = args
    idx += 1
  end
  return transforms
end

@unroll function __merge_transforms(transforms, i, coords, cur_s, cur_t_ref)                                                      
  @unroll for transform in transforms                                                                                            
    (transform.kernel)(i, coords, cur_s, cur_t_ref, transform.args...)                                                           
  end                                                                                                                            
  return nothing                                                                                                               
end          

@unroll function _execute_callbacks_with_point(i, callbacks, coords, p)
  @unroll for callback in callbacks
    callback(i, coords, p)
  end
  return nothing
end

execute_callbacks(i, coords, cur_s, cur_t_ref, cur_step) = _execute_callbacks(i, coords, coords.callbacks, cur_s, cur_t_ref, cur_step)

@unroll function _execute_callbacks(i, coords, callbacks, cur_s, cur_t_ref, cur_step)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, cur_t_ref, cur_step)
  end
  return nothing
end