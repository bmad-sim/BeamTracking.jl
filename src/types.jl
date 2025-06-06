"""
    MemoryLayout

Abstract type for memory layout strategies. Two implementations are provided:
- `AoS`: Array of Structures
- `SoA`: Structure of Arrays
"""
abstract type MemoryLayout end
struct AoS <: MemoryLayout end
struct SoA <: MemoryLayout end

"""
    Bunch{A<:MemoryLayout,S,T}

Structure representing a particle bunch.

# Fields
- `species::Species`: Particle species (e.g., ELECTRON, PROTON)
- `Brho_ref::S`: Reference magnetic rigidity
- `v::T`: Matrix of particle coordinates
          First index is particle, second is coordinate (x, px, y, py, z, pz)
          px, py are normalized momenta, pz is momentum deviation
"""
mutable struct Bunch{A<:MemoryLayout,S,T}
  species::Species   # Species
  Brho_ref::S        # Reference magnetic rigidity, used fornormalization of phase space coordinates
  const v::T         # Matrix of particle coordinates
  function Bunch{A}(species, Brho_ref, v) where {A}
    return new{A,typeof(Brho_ref),typeof(v)}(species, Brho_ref, v)
  end
end

# Constants for coordinate indexing
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6

"""
    soaview(bunch::Bunch{A}) where {A}

Get a Structure of Arrays view of the particle coordinates.
"""
soaview(bunch::Bunch{A}) where {A} = A == AoS ? transpose(bunch.v) : bunch.v

"""
    aosview(bunch::Bunch{A}) where {A}

Get an Array of Structures view of the particle coordinates.
"""
aosview(bunch::Bunch{A}) where {A} = A == AoS ? bunch.v : transpose(bunch.v)

"""
    get_N_particle(bunch::Bunch{A}) where {A}

Get the number of particles in the bunch.
"""
get_N_particle(bunch::Bunch{A}) where {A} = A == AoS ? size(bunch.v, 2) : size(bunch.v, 1)

"""
    setproperty!(bunch::Bunch{A,S}, key::Symbol, value) where {A,S}

Update bunch properties, handling special cases for Brho_ref and species changes.
Automatically adjusts particle momenta when Brho_ref or species changes.

# Arguments
- `bunch`: The particle bunch to modify
- `key`: Property to update (:Brho_ref or :species)
- `value`: New value for the property
"""
function setproperty!(bunch::Bunch{A,S}, key::Symbol, value) where {A,S}
  if key == :Brho_ref
    if value == bunch.Brho_ref
      return value
    end
    v = soaview(bunch)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, value)
    setfield!(bunch, :Brho_ref, S(value))
  elseif key == :species
    if value == bunch.species
      return value
    end
    v = soaview(bunch)
    Brho_final = bunch.Brho_ref*chargeof(bunch.species)/chargeof(value)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, Brho_final)
    setfield!(bunch, :Brho_ref, S(Brho_final))
    setfield!(bunch, :species, value)
  else
    setfield!(bunch, key, value)
  end
end

"""
    Bunch(N::Integer; mem=SoA, Brho_ref=NaN, species=ELECTRON)

Create a new bunch with N particles.

# Arguments
- `N`: Number of particles
- `mem`: Memory layout (SoA or AoS)
- `Brho_ref`: Reference magnetic rigidity
- `species`: Particle species

# Returns
A new `Bunch` instance with randomly initialized coordinates
"""
function Bunch(N::Integer; mem=SoA, Brho_ref=NaN, species=ELECTRON)
  if mem == SoA
    return Bunch{mem}(species, Brho_ref, rand(N,6))
  elseif mem == AoS
    return Bunch{mem}(species, Brho_ref, rand(6,N))
  else
    error("Invalid memory layout specification")
  end
end

"""
    Bunch(v::AbstractArray; mem=SoA, Brho_ref=NaN, species=ELECTRON)

Create a new bunch from existing coordinates.

# Arguments
- `v`: Matrix of particle coordinates
- `mem`: Memory layout (SoA or AoS)
- `Brho_ref`: Reference magnetic rigidity
- `species`: Particle species

# Returns
A new `Bunch` instance with the provided coordinates
"""
function Bunch(v::AbstractArray; mem=SoA, Brho_ref=NaN, species=ELECTRON)
  if mem == SoA
    size(v, 2) == 6 || error("For SoA the number of columns must be equal to 6")
  elseif mem == AoS
    size(v, 1) == 6 || error("For SoA the number of rows must be equal to 6")
  else
    error("Invalid memory layout specification")
  end
  return Bunch{mem}(species, Brho_ref, v)
end

"""
    ParticleView{S,T}

View into a single particle within a bunch.

# Fields
- `species::Species`: Particle species
- `Brho_ref::S`: Reference magnetic rigidity
- `index::Int`: Particle index
- `v::T`: View of particle coordinates
"""
struct ParticleView{S,T}
  species::Species
  Brho_ref::S     
  index::Int
  v::T    
end

"""
    ParticleView(bunch::Bunch{A}, i=1) where {A}

Create a view of a single particle in the bunch.

# Arguments
- `bunch`: The particle bunch
- `i`: Index of the particle to view (default: 1)

# Returns
A `ParticleView` instance for the specified particle
"""
function ParticleView(bunch::Bunch{A}, i=1) where {A}
  v = aosview(bunch)
  return ParticleView(bunch.species, bunch.Brho_ref, i, view(v, :, i))
end