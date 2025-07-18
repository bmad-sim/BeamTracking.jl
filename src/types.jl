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

struct Bunch{B,S,V,Q}
  species::Species # Species
  Brho_ref::B      # 0D Array: Defines normalization of phase space coordinates
  state::S   # Array of particle states
  v::V       # Matrix of particle coordinates
  q::Q       # Matrix of particle quaternions if spin else nothing 
  function Bunch(species, Brho_ref, state, v, q=nothing)
    return new{typeof(Brho_ref),typeof(state),typeof(v),typeof(q)}(species, Brho_ref, state, v, q)
  end
end

# Necessary for GPU compatibility:
Adapt.@adapt_structure Bunch


get_N_particle(bunch::Bunch) = size(bunch.v, 1)

function Bunch(N::Integer; Brho_ref=NaN, species=ELECTRON, spin=false)
  v = rand(N,6)
  q = spin ? rand(N,4) : nothing
  state = similar(v, State.T, N)
  state .= State.Alive
  return Bunch(species, fill(Brho_ref), state, v, q)
end

function Bunch(v::AbstractArray, q=nothing; Brho_ref=NaN, species=ELECTRON)
  size(v, 2) == 6 || error("The number of columns must be equal to 6")
  N_particle = size(v, 1)
  state = similar(v, State.T, N_particle)
  state .= State.Alive
  return Bunch(species, fill(Brho_ref), state, v, q)
end

# Update momenta for change to Brho_ref or change to species
function setproperty!(bunch::Bunch{B}, key::Symbol, value) where {B}
  if key == :Brho_ref
    error("Updating reference energy of bunch calculation not yet implemented")
    #=
    if value == bunch.Brho_ref
      return value
    end
    v = soaview(bunch)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, value)
    setfield!(bunch, :Brho_ref, B(value))
    =#
  elseif key == :species
    error("Updating species of bunch (which affects Brho_ref) not yet implemented")
    #=
    if value == bunch.species
      return value
    end
    v = soaview(bunch)
    Brho_final = bunch.Brho_ref*chargeof(bunch.species)/chargeof(value)
    launch!(Exact.update_P0!, v, nothing, bunch.Brho_ref, Brho_final)
    setfield!(bunch, :Brho_ref, B(Brho_final))
    setfield!(bunch, :species, value)
    =#
  else
    setfield!(bunch, key, value)
  end
end