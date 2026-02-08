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
struct Coords{S,V,Q}
  state::S # Array of particle states
  v::V     # Matrix of particle coordinates
  q::Q     # Matrix of particle quaternions if spin else nothing 
  function Coords(state, v, q)
    if !isnothing(q) && eltype(v) != eltype(q)
      error("Cannot initialize Coords with orbital coordinates of type $(eltype(v))
             and quaternion coordinates of type $(typeof(q)).")
    end
    return new{typeof(state),typeof(v),typeof(q)}(state, v, q)
  end
end

mutable struct Bunch{B,T,C<:Coords}
  species::Species # Species
  p_over_q_ref::B         # Defines normalization of phase space coordinates
  t_ref::T         # Reference time
  const coords::C  # GPU compatible structure of particles
end

# Necessary for GPU compatibility:
Adapt.@adapt_structure Coords

get_N_particle(bunch::Bunch) = size(bunch.coords.v, 1)

function Bunch(N::Integer; p_over_q_ref=NaN, t_ref=0., species=Species(), spin=false)
  v = rand(N,6)
  q = spin ? rand(N,4) : nothing
  state = similar(v, UInt8, N)
  state .= STATE_ALIVE
  return Bunch(species, p_over_q_ref, t_ref, Coords(state, v, q))
end

function Bunch(v::AbstractMatrix, q=nothing; p_over_q_ref=NaN, t_ref=0., species=Species())
  size(v, 2) == 6 || error("The number of columns must be equal to 6")
  N_particle = size(v, 1)
  state = similar(v, UInt8, N_particle)
  state .= STATE_ALIVE
  return Bunch(species, p_over_q_ref, t_ref, Coords(state, v, q))
end

function Bunch(v::AbstractVector, q=nothing; p_over_q_ref=NaN, t_ref=0., species=Species())
  length(v) == 6 || error("Bunch accepts a N x 6 matrix of N particle coordinates,
                            or alternatively a single particle as a vector. Received 
                            a vector of length $(length(v))")
  return Bunch(reshape(v, (1,6)), q; p_over_q_ref=p_over_q_ref, t_ref=t_ref, species=species)
end

struct ParticleView{B,T,S,V,Q}
  index::Int
  species::Species
  p_over_q_ref::B  
  t_ref::T   
  state::S
  v::V
  q::Q    
  ParticleView(args...) = new{typeof.(args)...}(args...)
end

function ParticleView(bunch::Bunch, i=1)
  v = bunch.coords.v
  q = bunch.coords.q
  return ParticleView(i, bunch.species, bunch.p_over_q_ref, bunch.t_ref, bunch.coords.state[i], view(v, :, i), isnothing(q) ? q : view(q, :, i))
end