#=
Example callback:

function mycallback(i, coords, cur_s, cur_t_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
  # Transform out of body frame:
  transforms_out!(i, coords, cur_s, cur_t_ref)

  # Do stuff in global frame...

  # Transform back into body frame:
  transforms_in!(i, coords, cur_s, cur_t_ref)
end

=#

# Constructs a single callback with signature (i, coords, cur_s, cur_t_ref) that executes all 
# callbacks in the callbacks tuple, including the transforms_out and transforms_in evaluated at 
# the given reference time and energy.
function construct_main_callback(coords, _transforms_out, _transforms_in, t_ref_transform, beta_gamma_ref_transform, ds_step, g)
  # Note: in keeping consistency with the rest of TimeDependentParam, transforms with time-dependent params 
  # (e.g. misalignment with time-dependent x_offset) only evaluated at entrance of element. This is actually 
  # needed to ensure everything is consistent, i.e., transform at exit consistent with rest of bunch. 
  # perhaps we will want to refactor this if we want TimeDependentParam to evolve inside of elements, 
  # but much lower priority.
  let _transforms_in = _transforms_in, _transforms_out = _transforms_out, t_ref_transform=t_ref_transform, beta_gamma_ref_transform=beta_gamma_ref_transform, callbacks=coords.callbacks, ds_step=ds_step, g=g
    return (i, coords, cur_s, cur_dt_ref) -> begin
      transforms_out! = _merge_transforms(_evaluate_transforms_args(i, coords, _transforms_out, t_ref_transform, beta_gamma_ref_transform))
      transforms_in! = _merge_transforms(_evaluate_transforms_args(i, coords, _transforms_in, t_ref_transform, beta_gamma_ref_transform))
      _execute_callbacks(i, callbacks, coords, cur_s, t_ref_transform+cur_dt_ref, ds_step, g, transforms_out!, transforms_in!)
      return nothing
    end
  end
end

@unroll function _evaluate_transforms_args(i, coords, transforms, t_ref_transform, beta_gamma_ref_transform)
  idx = 1
  @unroll for transform in transforms
    bargs = process_batch_args(i, transform.args)
    args = process_time_args(i, coords, bargs, t_ref_transform, beta_gamma_ref_transform)
    @reset transforms[idx].args = args
    idx += 1
  end
  return transforms
end

function _merge_transforms(transforms)                                                                                           
  let transforms=transforms                                                                                                      
    return (i, coords, cur_s, cur_t_ref) -> __merge_transforms(transforms, i, coords, cur_s, cur_t_ref)                           
  end                                                                                                                            
end                                                                                                                            
                                                                                                                                 
@unroll function __merge_transforms(transforms, i, coords, cur_s, cur_t_ref)                                                      
  @unroll for transform in transforms                                                                                            
    (transform.kernel)(i, coords, cur_s, cur_t_ref, transform.args...)                                                           
  end                                                                                                                            
  return nothing                                                                                                               
end          

@unroll function _execute_callbacks(i, callbacks, coords, cur_s, cur_t_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
  @unroll for callback in callbacks
    callback(i, coords, cur_s, cur_t_ref, last_ds_step, last_g, transforms_out!, transforms_in!)
  end
  return nothing
end