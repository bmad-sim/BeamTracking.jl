# Hook for effects that act on the whole bunch between slices of the
# elements of a Beamline (e.g. a wake computed from all particles).
# Used through `track!(bunch, bl; effects = (effect, ...))`.

"""
    SlicePoint{T,S}

Entrance of one slice of an element, where the `effects` of `track!`
are applied (see `effect_kick!`).

# Fields
- `ix`: index of the element in the `line` of the `Beamline`.
- `k`: slice number, `1` to `n_slices`.
- `n_slices`: number of slices of the element.
- `s0`: slice entrance from the element entrance [m].
- `L_slice`: slice length [m].
- `L`: element length [m].
- `s`: slice entrance from the `Beamline` entrance [m] (`bunch.s`
  there, of the type of `bunch.t_ref`).
"""
struct SlicePoint{T,S}
  ix::Int
  k::Int
  n_slices::Int
  s0::T
  L_slice::T
  L::T
  s::S
end

"""
    effect_setup(effect, bunch, bl; kwargs...) -> state

Called once per `track!(bunch, bl; effects)` for each effect, before
the first element is tracked, with the bunch's species and reference
already resolved. The returned `state` is passed to every
`effect_kick!` call of that effect. Packages add methods for their
own effect types; there is no default method.

The keywords are those of the `track!` call except `notify` and
`effects` (`scalar_params`, `ramp_particle_energy_without_rf`,
`ramp_update_each_particle`, `rf_on`, `batch_start` and any other
keyword given to it), so that an effect that tracks particles itself
(e.g. a reference orbit) can track them as the call does. Methods
must accept them (`; kwargs...`).
"""
function effect_setup end

"""
    effect_kick!(effect, state, bunch, ele, sp::SlicePoint)

Called at the entrance of every slice of every element (an element
with `n_slices = 1` is one slice), before the slice is tracked. The
method modifies `bunch` in place: its coordinates are at the slice
entrance, and `bunch.p_over_q_ref` and `bunch.t_ref` are those left
by the previous slice. Effects are applied in the order of the
`effects` tuple. Packages add methods for their own effect types;
there is no default method.

# Arguments
- `effect`: the entry of `effects`.
- `state`: what `effect_setup` returned for it.
- `bunch`: the tracked `Bunch`.
- `ele`: the element being tracked.
- `sp`: the slice (`SlicePoint`).
"""
function effect_kick! end

# Applies each effect with its state, in tuple order
function apply_effects!(effects::Tuple, states::Tuple, bunch, ele, sp)
  foreach(effects, states) do effect, state
    effect_kick!(effect, state, bunch, ele, sp)
  end
  return nothing
end
