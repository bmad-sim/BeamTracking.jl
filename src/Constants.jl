module Constants
using AtomicAndPhysicalConstants
@APCdef tupleflag=false
isnullspecies(species::Species) = getfield(species, :kind) == Kind.NULL
end