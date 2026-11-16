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

- `s`: s from the beginning of the `Beamline` (from the beginning of the element when tracking a
  single element)
- `t_ref`: reference time
- `step`: integration steps from the beginning of the `Beamline` (`n_slices*n_steps` per element,
  the last one at the exit of the element; from the beginning of the element when tracking a single
  element)
- `beta_gamma_ref`: reference energy (beta*gamma)
- `last_ds_step`: length of the last integration step
- `last_g`: curvature `(gx, gy)` of the last integration step

`transforms_out!(i, coords, p)` and `transforms_in!(i, coords, p)` transform particle `i` out of 
and back into the element body frame (e.g. of a misaligned element).

`t_ref` and the transforms are only computed when used.
"""
struct TrackPoint{S,SE,T,D,B,L,G,TOUT,TIN}
  s::S                # from the Beamline entrance
  s_ele::SE           # from the element entrance (for the transforms)
  t_enter::T          # reference time at the launch entrance
  dt::D               # reference time since then, a number or a _DtRef (computed on access)
  step::Int
  beta_gamma_ref::B
  last_ds_step::L
  last_g::G
  batch_start::Int
  transforms_out::TOUT
  transforms_in::TIN
end

# The reference time from the launch entrance to s in the integration of kernel ker with params,
# computed only when the time is used (it costs a square root and divisions at every step)
struct _DtRef{S,K,P}
  s::S
  ker::K
  params::P
end

_dt_value(dt) = dt
_dt_value(dt::_DtRef) = compute_dt_ref(dt.s, dt.ker, dt.params)

Base.getproperty(p::TrackPoint, key::Symbol) =
  key === :t_ref ? getfield(p, :t_enter) + _dt_value(getfield(p, :dt)) : getfield(p, key)

Base.propertynames(::TrackPoint) = (:s, :t_ref, :step, :beta_gamma_ref, :last_ds_step, :last_g)

# Note: in keeping consistency with the rest of TimeDependentParam, transforms with time-dependent params 
# (e.g. misalignment with time-dependent x_offset) only evaluated at entrance of element. This is actually 
# needed to ensure everything is consistent, i.e., transform at exit consistent with rest of bunch. 
# perhaps we will want to refactor this if we want TimeDependentParam to evolve inside of elements, 
# but much lower priority.
transforms_out!(i, coords, p::TrackPoint) = _transform!(i, coords, p, getfield(p, :transforms_out))
transforms_in!(i, coords, p::TrackPoint) = _transform!(i, coords, p, getfield(p, :transforms_in))

"""
    has_transforms(p::TrackPoint)

`true` if `transforms_out!` and `transforms_in!` change the coordinates at `p` (e.g. in a 
misaligned element), `false` if they do nothing. Known at compile time, so that code for the 
transforms can be left out where there are none.
"""
@generated function has_transforms(::TrackPoint{S,SE,T,D,B,L,G,TOUT,TIN}) where {S,SE,T,D,B,L,G,TOUT,TIN}
  blank(K) = K <: KernelCall{typeof(blank_kernel!)}
  return !(all(blank, fieldtypes(TOUT)) && all(blank, fieldtypes(TIN)))
end

function _transform!(i, coords, p, transforms)
  t_enter = getfield(p, :t_enter)
  args = _evaluate_transforms_args(i, coords, getfield(p, :batch_start), transforms, t_enter,
                                   getfield(p, :beta_gamma_ref))
  __merge_transforms(args, i, coords, getfield(p, :s_ele), p.t_ref)
  return nothing
end

"""
    runs_at(callback, step) -> Bool

`false` if `callback` does nothing at the integration step `step` (counted as `p.step` of the 
`TrackPoint`), so that it does not need to be called there. Defaults to `true`. A callback that runs 
only at some steps can define a method: the integration then only stops at the steps where it runs 
(and elements where it does not run at the exit skip it altogether), instead of keeping the code to 
call it in every integration step. It must be cheap, and the same for all particles.
"""
runs_at(callback, step) = true
@inline runs_at(callbacks::Tuple{}, step) = false
@inline runs_at(callbacks::Tuple, step) =
  runs_at(first(callbacks), step) || runs_at(Base.tail(callbacks), step)

# A single callback with signature (i, coords, cur_s, cur_dt_ref, cur_step) that executes all
# callbacks in the callbacks tuple with the TrackPoint. cur_s and cur_step are from the launch 
# entrance: the callbacks get them from the Beamline entrance (s_enter, the s of the launch 
# entrance, and step0, the steps before it, are added), and the transforms get s from the element 
# entrance (s0, the s of the launch entrance within the element, nonzero for a slice of an element, 
# is added). cur_dt_ref is the reference time from the launch entrance (or a _DtRef).
struct MainCallback{C,TOUT,TIN,T,B,D,G,S0,SE}
  callbacks::C
  batch_start::Int
  transforms_out::TOUT
  transforms_in::TIN
  t_ref_transform::T
  beta_gamma_ref_transform::B
  ds_step::D
  g::G
  s0::S0
  s_enter::SE
  step0::Int
end

function (m::MainCallback)(i, coords, cur_s, cur_dt_ref, cur_step)
  p = TrackPoint(m.s_enter+cur_s, m.s0+cur_s, m.t_ref_transform, cur_dt_ref, m.step0+cur_step, m.beta_gamma_ref_transform, m.ds_step, m.g, m.batch_start, m.transforms_out, m.transforms_in)
  _execute_callbacks_with_point(i, m.callbacks, coords, p)
  return nothing
end

@inline runs_at(m::MainCallback, cur_step) = runs_at(m.callbacks, m.step0 + cur_step)

construct_main_callback(coords, batch_start, _transforms_out, _transforms_in, t_ref_transform, beta_gamma_ref_transform, ds_step, g, s0=0, s_enter=s0, step0=0) =
  (MainCallback(coords.callbacks, batch_start, _transforms_out, _transforms_in, t_ref_transform, beta_gamma_ref_transform, ds_step, g, s0, s_enter, step0),)

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