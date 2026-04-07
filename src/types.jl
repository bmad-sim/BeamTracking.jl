const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6
const Q0  = 1
const QX  = 2
const QY  = 3
const QZ  = 4

const STATE_PREBORN    = UInt8(0)
const STATE_ALIVE      = UInt8(1)
const STATE_LOST       = UInt8(2)
const STATE_LOST_NEG_X = UInt8(3)
const STATE_LOST_POS_X = UInt8(4)
const STATE_LOST_NEG_Y = UInt8(5)
const STATE_LOST_POS_Y = UInt8(6)
const STATE_LOST_PZ    = UInt8(7)
const STATE_LOST_Z     = UInt8(8)

# Always SOA
struct Coords{S,V,Q,W,T}
  state::S # Array of particle states
  v::V     # Matrix of particle coordinates
  q::Q     # Matrix of particle quaternions if spin else nothing 
  weight::W     # Array of particle weights if weighted else nothing
  callbacks::T  # Tuple of functions to evaluate inside kernels
  function Coords(state, v, q, weight, callbacks)
    if !isnothing(q) && eltype(v) != eltype(q)
      error("Cannot initialize Coords with orbital coordinates of type $(eltype(v))
             and quaternion coordinates of type $(typeof(q)).")
    end
    return new{typeof(state),typeof(v),typeof(q),typeof(weight),typeof(callbacks)}(state, v, q, weight, callbacks)
  end
end

mutable struct Bunch{B,T,C<:Coords}
  species::Species # Species
  p_over_q_ref::B         # Defines normalization of phase space coordinates
  t_ref::T         # Reference time
  const coords::C  # GPU compatible structure of particles
end

function Base.getproperty(b0::Bunch, key::Symbol)
  if key in (:state, :v, :q, :weight, :callbacks)
    return getproperty(b0.coords, key)
  else
    return getfield(b0, key)
  end
end

Base.propertynames(b0::Bunch) = (:state, :v, :q, :weight, :callbacks, :species, :p_over_q_ref, :t_ref, :coords)

# Necessary for GPU compatibility:
Adapt.@adapt_structure Coords

get_N_particle(bunch::Bunch) = size(bunch.coords.v, 1)

"""
    Bunch(; v, state, spin, q, weight, callbacks, p_over_q_ref, t_ref, species)

Construct a `Bunch` of particles for tracking.

# Required Keyword Arguments
- `v::AbstractMatrix`: Matrix of particle phase space coordinates with shape 
  `(N, 6)`, where `N` is the number of particles and the 6 columns are the 
  canonical coordinates `[x px y py z pz]`.

# Optional Keyword Arguments
- `state`: Array of particle states with length `N`. Each element is a `UInt8` 
  flag indicating the particle status (e.g. `0x1 == STATE_ALIVE`). Defaults to all 
  particles alive.
- `spin::Bool=false`: If `true`, allocate identity quaternions for spin tracking.
- `q`: Matrix of spin quaternions with shape `(N, 4)`, where each row is a 
  unit quaternion `[q0 q1 q2 q3]`. Element type must match that of `v`.
- `weight`: Array of particle weights with length `N`, or `nothing` if 
   even weights for all particles. Defaults to `nothing`.
- `callbacks`: Tuple of functions to evaluate inside tracking kernels, which should 
   be of the form `callback(coords::Coords, ds_step, g)` where `coords` is the structure 
   storing particle coordinates AFTER taking an integration step of size `ds_step`, and 
   `g` is a tuple of the instantaneous coordinate system curvatures in x and y, i.e. 
   `g = (gx, gy)`. Note that if more information is needed (e.g. previous step `coords` 
   or perhaps `s` position), then the callback should be a *closure* that reads/writes 
   to an outside mutable state each execution. 
  Defaults to `()`.
- `p_over_q_ref`: Reference momentum-over-charge defining the normalization 
  of the phase space coordinates.
- `t_ref=0.`: Reference time.
- `species=Species()`: Particle species.

# Example
```julia
N = 1000
v = zeros(N, 6)
bunch = Bunch(v=v, species=Species("electron"), p_over_q_ref=-60.0)

# With spin tracking:
bunch = Bunch(v=v, spin=true, species=Species("electron"), p_over_q_ref=-60.0)
```
"""
function Bunch(;
  v::AbstractMatrix,
  state=(s = similar(v, UInt8, size(v, 1)); s .= STATE_ALIVE; s),
  spin::Bool=false,
  q= spin ? (qs = similar(v, (size(v, 1), 4)); qs .= 0; qs[:,1] .= 1; qs) : nothing,
  weight=nothing,
  callbacks=(),
  p_over_q_ref=NaN, 
  t_ref=0., 
  species=Species(),
)
  size(v, 2) == 6 || error("The number of columns of the particle coordinates vector `v` must be equal to 6")
  return Bunch(species, p_over_q_ref, t_ref, Coords(state, v, q, weight, callbacks))
end

function Bunch(v::AbstractMatrix, q=nothing, weight=nothing; p_over_q_ref=NaN, t_ref=0., species=Species(), callbacks=())
  size(v, 2) == 6 || error("The number of columns must be equal to 6")
  N_particle = size(v, 1)
  state = similar(v, UInt8, N_particle)
  state .= STATE_ALIVE
  return Bunch(species, p_over_q_ref, t_ref, Coords(state, v, q, weight, callbacks))
end

function Bunch(v::AbstractVector, q=nothing, weight=nothing; p_over_q_ref=NaN, t_ref=0., species=Species(), callbacks=())
  length(v) == 6 || error("Bunch accepts a N x 6 matrix of N particle coordinates,
                            or alternatively a single particle as a vector. Received 
                            a vector of length $(length(v))")
  return Bunch(reshape(v, (1,6)), q, weight; p_over_q_ref=p_over_q_ref, t_ref=t_ref, species=species, callbacks=callbacks)
end

struct ParticleView{B,T,S,V,Q,W}
  index::Int
  species::Species
  p_over_q_ref::B  
  t_ref::T   
  state::S
  v::V
  q::Q
  weight::W
  ParticleView(args...) = new{typeof.(args)...}(args...)
end

function ParticleView(bunch::Bunch, i=1)
  v = bunch.coords.v
  q = bunch.coords.q
  weight = bunch.coords.weight
  return ParticleView(i, bunch.species, bunch.p_over_q_ref, bunch.t_ref, bunch.coords.state[i], view(v, :, i), isnothing(q) ? q : view(q, :, i), isnothing(weight) ? weight : weight[i])
end