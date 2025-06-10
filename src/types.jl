abstract type MemoryLayout end
struct AoS <: MemoryLayout end
struct SoA <: MemoryLayout end

@enumx State::UInt8 Preborn Alive Lost Lost_Neg_X Lost_Pos_X Lost_Neg_Y Lost_Pos_Y Lost_Pz Lost_Z     

mutable struct Bunch{mem<:MemoryLayout,B,S,V,Q}
  species::Species # Species
  Brho_ref::B      # Defines normalization of phase space coordinates
  const state::S   # Array of particle states
  const v::V       # Matrix of particle coordinates
  const q::Q       # Matrix of particle quaternions if spin else nothing 
  function Bunch{mem}(species, Brho_ref, state, v, q=nothing) where {mem}
    return new{mem,typeof(Brho_ref),typeof(state),typeof(v),typeof(q)}(species, Brho_ref, state, v, q)
  end
end

unpack(bunch::Bunch) = bunch.state, soaview(bunch)...

# Index particle i coordinate x as (i,1) , px as (i,2), etc
function soaview(bunch::Bunch{mem}) where {mem}
  if mem == AoS
    v = transpose(bunch.v)
    q = isnothing(bunch.q) ? nothing : transpose(bunch.q)
  else
    v = bunch.v
    q = isnothing(bunch.q) ? nothing : bunch.q
  end
  return v, q
end

function aosview(bunch::Bunch{mem}) where {mem}
  if mem == AoS
    v = bunch.v
    q = Q == Nothing ? nothing : bunch.q
  else
    v = transpose(bunch.v)
    q = Q == Nothing ? nothing : transpose(bunch.q)
  end
  return v, q
end

get_N_particle(bunch::Bunch{mem}) where {mem} = mem == AoS ? size(bunch.v, 2) : size(bunch.v, 1)

function Bunch(N::Integer; mem=SoA, Brho_ref=NaN, species=ELECTRON, spin=false)
  if mem == SoA
    v = rand(N,6)
    q = spin ? rand(N,4) : nothing
  else
    v = rand(6,N)
    q = spin ? rand(4,N) : nothing
  end
  state = similar(v, State.T, N)
  state .= State.Preborn
  return Bunch{mem}(species, Brho_ref, state, v, q)
end

function Bunch(v::AbstractArray, q=nothing; mem=SoA, Brho_ref=NaN, species=ELECTRON)
  if mem == SoA
    size(v, 2) == 6 || error("For SoA the number of columns must be equal to 6")
    N_particle = size(v, 1)
  else
    size(v, 1) == 6 || error("For SoA the number of rows must be equal to 6")
    N_particle = size(v, 2)
  end
  state = similar(v, State.T, N_particle)
  state .= State.Preborn
  return Bunch{mem}(species, Brho_ref, state, v, q)
end

struct ParticleView{mem<:MemoryLayout,B,S,V,Q}
  index::Int
  species::Species
  Brho_ref::B     
  state::S
  v::V
  q::Q    
end

function ParticleView(bunch::Bunch, i=1)
  v, q = aosview(bunch)
  return ParticleView(i, bunch.species, bunch.Brho_ref, bunch.state[i], view(v, :, i), isnothing(q) ? q : view(q, :, i))
end

# Update momenta for change to Brho_ref or change to species
function setproperty!(bunch::Bunch{mem,B}, key::Symbol, value) where {mem,B}
  if key == :Brho_ref
    error("Updating reference energy of bunch calculation not yet implemented")
    if value == bunch.Brho_ref
      return value
    end
    v,__ = soaview(bunch)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, value)
    setfield!(bunch, :Brho_ref, B(value))
  elseif key == :species
    error("Updating species of bunch (which affect Brho_ref) not yet implemented")
    if value == bunch.species
      return value
    end
    v,__ = soaview(bunch)
    Brho_final = bunch.Brho_ref*chargeof(bunch.species)/chargeof(value)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, Brho_final)
    setfield!(bunch, :Brho_ref, B(Brho_final))
    setfield!(bunch, :species, value)
  else
    setfield!(bunch, key, value)
  end
end
