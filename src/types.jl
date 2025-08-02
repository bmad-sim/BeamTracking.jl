abstract type MemoryLayout end
struct AoS <: MemoryLayout end
struct SoA <: MemoryLayout end
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6

@enumx State::UInt8 Preborn Alive Lost Lost_Neg_X Lost_Pos_X Lost_Neg_Y Lost_Pos_Y Lost_Pz Lost_Z   

mutable struct Bunch{mem<:MemoryLayout,B,S,V,Q}
  species::Species # Species
  rigidity::B      # Defines normalization of phase space coordinates
  const state::S   # Array of particle states
  const v::V       # Matrix of particle coordinates
  const q::Q       # Matrix of particle quaternions if spin else nothing 
  function Bunch{mem}(species, rigidity, state, v, q=nothing) where {mem}
    return new{mem,typeof(rigidity),typeof(state),typeof(v),typeof(q)}(species, rigidity, state, v, q)
  end
end

# Always SOA
struct BunchView{S,V,Q}
  state::S
  v::V
  q::Q
  function BunchView(b::Bunch)
    state = b.state
    v = soaview(b)
    q = soaviewq(b)
    return new{typeof.((state,v,q))...}(state,v,q)
  end
  BunchView(args...) = new{typeof.(args)...}(args...)
end

# Necessary for GPU compatibility:
Adapt.@adapt_structure BunchView

# Index particle i coordinate x as (i,1) , px as (i,2), etc
soaview(bunch::Bunch{A}) where {A} = A == AoS ? transpose(bunch.v) : bunch.v
aosview(bunch::Bunch{A}) where {A} = A == AoS ? bunch.v : transpose(bunch.v)

soaviewq(bunch::Bunch{A}) where {A} = (A == SoA || isnothing(bunch.q)) ? bunch.q : transpose(bunch.q)
aosviewq(bunch::Bunch{A}) where {A} = (A == AoS || isnothing(bunch.q)) ? bunch.q : transpose(bunch.q)

get_N_particle(bunch::Bunch{A}) where {A} = A == AoS ? size(bunch.v, 2) : size(bunch.v, 1)

function Bunch(N::Integer; mem=SoA, rigidity=NaN, species=ELECTRON, spin=false)
  if mem == SoA
    v = rand(N,6)
    q = spin ? rand(N,4) : nothing
  else
    v = rand(6,N)
    q = spin ? rand(4,N) : nothing
  end
  state = similar(v, State.T, N)
  state .= State.Alive
  return Bunch{mem}(species, rigidity, state, v, q)
end

function Bunch(v::AbstractArray, q=nothing; mem=SoA, rigidity=NaN, species=ELECTRON)
  if mem == SoA
    size(v, 2) == 6 || error("For SoA the number of columns must be equal to 6")
    N_particle = size(v, 1)
  else
    size(v, 1) == 6 || error("For SoA the number of rows must be equal to 6")
    N_particle = size(v, 2)
  end
  state = similar(v, State.T, N_particle)
  state .= State.Alive
  return Bunch{mem}(species, rigidity, state, v, q)
end

struct ParticleView{B,S,V,Q}
  index::Int
  species::Species
  rigidity::B     
  state::S
  v::V
  q::Q    
  ParticleView(args...) = new{typeof.(args)...}(args...)
end

function ParticleView(bunch::Bunch, i=1)
  v = aosview(bunch)
  q = aosviewq(bunch)
  return ParticleView(i, bunch.species, bunch.rigidity, bunch.state[i], view(v, :, i), isnothing(q) ? q : view(q, :, i))
end

# Update momenta for change to rigidity or change to species
function setproperty!(bunch::Bunch{mem,B}, key::Symbol, value) where {mem,B}
  if key == :rigidity
    error("Updating reference energy of bunch calculation not yet implemented")
    #=
    if value == bunch.rigidity
      return value
    end
    v = soaview(bunch)
    launch!(Exact.update_P0!, v, nothing, bunch.rigidity, value)
    setfield!(bunch, :rigidity, B(value))
    =#
  elseif key == :species
    error("Updating species of bunch (which affects rigidity) not yet implemented")
    #=
    if value == bunch.species
      return value
    end
    v = soaview(bunch)
    rigidity_final = bunch.rigidity*chargeof(bunch.species)/chargeof(value)
    launch!(Exact.update_P0!, v, nothing, bunch.rigidity, rigidity_final)
    setfield!(bunch, :rigidity, B(rigidity_final))
    setfield!(bunch, :species, value)
    =#
  else
    setfield!(bunch, key, value)
  end
end