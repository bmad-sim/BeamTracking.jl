abstract type MemoryLayout end
struct AoS <: MemoryLayout end
struct SoA <: MemoryLayout end

mutable struct Bunch{A<:MemoryLayout,S,T}
  species::Species # Species
  Brho_0::S        # Defines normalization of phase space coordinates
  const v::T       # Matrix of particle coordinates
  function Bunch{A}(species, Brho_0, v) where {A}
    return new{A,typeof(Brho_0),typeof(v)}(species, Brho_0, v)
  end
end

# Index particle i coordinate x as (i,1) , px as (i,2), etc
soaview(bunch::Bunch{A}) where {A} = A == AoS ? transpose(bunch.v) : bunch.v
aosview(bunch::Bunch{A}) where {A} = A == AoS ? bunch.v : transpose(bunch.v)
get_N_particle(bunch::Bunch{A}) where {A} = A == AoS ? size(bunch.v, 2) : size(bunch.v, 1)

# Update momenta for change to Brho_0 or change to species
function setproperty!(bunch::Bunch{A,S}, key::Symbol, value) where {A,S}
  if key == :Brho_0
    if value == bunch.Brho_0
      return value
    end
    v = soaview(bunch)
    launch!(General.update_P0!, v, nothing, bunch.Brho_0, value)
    setfield!(bunch, :Brho_0, S(value))
  elseif key == :species
    if value == bunch.species
      return value
    end
    v = soaview(bunch)
    Brho_final = bunch.Brho_0*chargeof(bunch.species)/chargeof(value)
    launch!(General.update_P0!, v, nothing, bunch.Brho_0, Brho_final)
    setfield!(bunch, :Brho_0, S(Brho_final))
    setfield!(bunch, :species, value)
  else
    setfield!(bunch, key, value)
  end
end

function Bunch(N::Integer; mem=SoA, Brho_0=60.0, species=ELECTRON)
  if mem == SoA
    return Bunch{mem}(species, Brho_0, rand(N,6))
  elseif mem == AoS
    return Bunch{mem}(species, Brho_0, rand(6,N))
  else
    error("Invalid memory layout specification")
  end
end

function Bunch(v::AbstractArray; mem=SoA, Brho_0=60.0, species=ELECTRON)
  if mem == SoA
    size(v, 2) == 6 || error("For SoA the number of columns must be equal to 6")
  elseif mem == AoS
    size(v, 1) == 6 || error("For SoA the number of rows must be equal to 6")
  else
    error("Invalid memory layout specification")
  end
  return Bunch{mem}(species, Brho_0, v)
end

struct ParticleView{S,T}
  species::Species
  Brho_0::S     
  index::Int
  v::T    
end

function ParticleView(bunch::Bunch{A}, i=1) where {A}
  v = aosview(bunch)
  return ParticleView(bunch.species, bunch.Brho_0, i, view(v, :, i))
end