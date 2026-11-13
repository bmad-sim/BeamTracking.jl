#=
Example callback:

function mycallback(i, coords, cur_s, cur_t_ref, cur_step, cur_beta_gamma_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
  # Transform out of body frame:
  transforms_out!(i, coords)

  # Do stuff in global frame...

  # Transform back into body frame:
  transforms_in!(i, coords)
end

=#

# Constructs a single callback with signature (i, coords, cur_s, cur_t_ref, cur_step) that executes 
# all callbacks in the callbacks tuple, including the transforms_out and transforms_in evaluated at 
# the given reference time and energy. cur_s and cur_step are from the launch entrance: the callbacks
# get them from the Beamline entrance (s_enter, the s of the launch entrance, and step0, the steps 
# before it, are added), and the transforms get s from the element entrance (s0, the s of the 
# launch entrance within the element, nonzero for a slice of an element, is added) and the time, 
# which they capture, so that they are called as transforms_out!(i, coords).
function construct_main_callback(coords, batch_start, _transforms_out, _transforms_in, t_ref_transform, beta_gamma_ref_transform, ds_step, g, s0=0, s_enter=s0, step0=0)
  # Note: in keeping consistency with the rest of TimeDependentParam, transforms with time-dependent params 
  # (e.g. misalignment with time-dependent x_offset) only evaluated at entrance of element. This is actually 
  # needed to ensure everything is consistent, i.e., transform at exit consistent with rest of bunch. 
  # perhaps we will want to refactor this if we want TimeDependentParam to evolve inside of elements, 
  # but much lower priority.
  let batch_start = batch_start, _transforms_in = _transforms_in, _transforms_out = _transforms_out, t_ref_transform=t_ref_transform, beta_gamma_ref_transform=beta_gamma_ref_transform, callbacks=coords.callbacks, ds_step=ds_step, g=g, s0=s0, s_enter=s_enter, step0=step0
    return ((i, coords, cur_s, cur_dt_ref, cur_step) -> begin
      cur_t_ref = t_ref_transform + cur_dt_ref
      transforms_out! = _merge_transforms(_evaluate_transforms_args(i, coords, batch_start, _transforms_out, t_ref_transform, beta_gamma_ref_transform), s0+cur_s, cur_t_ref)
      transforms_in! = _merge_transforms(_evaluate_transforms_args(i, coords, batch_start, _transforms_in, t_ref_transform, beta_gamma_ref_transform), s0+cur_s, cur_t_ref)
      _execute_callbacks_with_transforms(i, callbacks, coords, s_enter+cur_s, cur_t_ref, step0+cur_step, beta_gamma_ref_transform, ds_step, g, transforms_out!, transforms_in!)
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

function _merge_transforms(transforms, cur_s, cur_t_ref)
  let transforms=transforms, cur_s=cur_s, cur_t_ref=cur_t_ref
    return (i, coords) -> __merge_transforms(transforms, i, coords, cur_s, cur_t_ref)
  end
end                                                                                                                            
                                                                                                                                 
@unroll function __merge_transforms(transforms, i, coords, cur_s, cur_t_ref)                                                      
  @unroll for transform in transforms                                                                                            
    (transform.kernel)(i, coords, cur_s, cur_t_ref, transform.args...)                                                           
  end                                                                                                                            
  return nothing                                                                                                               
end          

@unroll function _execute_callbacks_with_transforms(i, callbacks, coords, cur_s, cur_t_ref, cur_step, cur_beta_gamma_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, cur_t_ref, cur_step, cur_beta_gamma_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
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