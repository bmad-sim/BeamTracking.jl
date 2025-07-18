using BeamTracking, Beamlines, GTPSA, BenchmarkTools

# Read in the Electron Storage Ring of the Electron-Ion Collider
include("test/lattices/esr.jl") # Beamline symbol is "ring"

# Currently only Linear tracking is supported:
foreach(t -> t.tracking_method = Linear(), ring.line)

# Construct a bunch:
N_particle = 100
b0 = Bunch(N_particle)

# Track the bunch through the ESR
track!(b0, ring)

# Also can track! individual elements
track!(b0, ring.line[1])

# And if you want to move the particle for-loop to the outside:
track!(b0, ring; outer_particle_loop=true)

# Can also track the bits representation:
bitsring = BitsBeamline(ring)
track!(b0, bitsring)

# GTPSA map:
const D = Descriptor(6, 1) # 6 variables to 1st order
v = @vars(D)
b0 = Bunch(v)

track!(b0, ring)
track!(b0, bitsring)