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

const State_Preborn    = UInt8(0)
const State_Alive      = UInt8(1)
const State_Lost       = UInt8(2)
const State_Lost_Neg_X = UInt8(3)
const State_Lost_Pos_X = UInt8(4)
const State_Lost_Neg_Y = UInt8(5)
const State_Lost_Pos_Y = UInt8(6)
const State_Lost_Pz    = UInt8(7)
const State_Lost_Z     = UInt8(8)

#@enumx State::UInt8 Preborn Alive Lost Lost_Neg_X Lost_Pos_X Lost_Neg_Y Lost_Pos_Y Lost_Pz Lost_Z   

# Always SOA
struct Coords{S,V,Q}
  state::S # Array of particle states
  v::V     # Matrix of particle coordinates
  q::Q     # Matrix of particle quaternions if spin else nothing 
end

mutable struct Bunch{B,C<:Coords}
  species::Species # Species
  R_ref::B         # Defines normalization of phase space coordinates
  const coords::C  # GPU compatible structure of particles
end

# Necessary for GPU compatibility:
Adapt.@adapt_structure Coords

get_N_particle(bunch::Bunch) = size(bunch.coords.v, 1)

function Bunch(N::Integer; R_ref=NaN, species=Species(), spin=false)
  v = rand(N,6)
  q = spin ? rand(N,4) : nothing
  state = similar(v, UInt8, N)
  state .= State_Alive
  return Bunch(species, R_ref, Coords(state, v, q))
end

function Bunch(v::AbstractMatrix, q=nothing; R_ref=NaN, species=Species())
  size(v, 2) == 6 || error("The number of columns must be equal to 6")
  N_particle = size(v, 1)
  state = similar(v, UInt8, N_particle)
  state .= State_Alive
  return Bunch(species, R_ref, Coords(state, v, q))
end

function Bunch(v::AbstractVector, q=nothing; R_ref=NaN, species=Species())
  length(v) == 6 || error("Bunch accepts a N x 6 matrix of N particle coordinates,
                            or alternatively a single particle as a vector. Received 
                            a vector of length $(length(v))")
  return Bunch(reshape(v, (1,6)), q; R_ref=R_ref, species=species)
end

struct ParticleView{B,S,V,Q}
  index::Int
  species::Species
  R_ref::B     
  state::S
  v::V
  q::Q    
  ParticleView(args...) = new{typeof.(args)...}(args...)
end

function ParticleView(bunch::Bunch, i=1)
  v = bunch.coords.v
  q = bunch.coords.q
  return ParticleView(i, bunch.species, bunch.R_ref, bunch.coords.state[i], view(v, :, i), isnothing(q) ? q : view(q, :, i))
end

# Update momenta for change to R_ref or change to species
function setproperty!(bunch::Bunch, key::Symbol, value)
  if key == :R_ref
    error("Updating reference energy of bunch calculation not yet implemented")
    #=
    if value == bunch.R_ref
      return value
    end
    v = soaview(bunch)
    launch!(Exact.update_P0!, v, nothing, bunch.R_ref, value)
    setfield!(bunch, :R_ref, B(value))
    =#
  elseif key == :species
    error("Updating species of bunch (which affects R_ref) not yet implemented")
    #=
    if value == bunch.species
      return value
    end
    v = soaview(bunch)
    R_ref_final = bunch.R_ref*chargeof(bunch.species)/chargeof(value)
    launch!(Exact.update_P0!, v, nothing, bunch.R_ref, R_ref_final)
    setfield!(bunch, :R_ref, B(R_ref_final))
    setfield!(bunch, :species, value)
    =#
  else
    setfield!(bunch, key, value)
  end
end