using CUDA, BeamTracking, Beamlines, PhysicalConstants

import  PhysicalConstants.CODATA2022: c_0 as c, m_e as m

include("simple_ring.jl") # Beamline symbol is "ring"
# Currently only Linear tracking is supported, enable it for each element
foreach(t -> t.tracking_method = Linear(), ring.line)
bitsring = BitsBeamline(ring)

# n_particles = parse(Int, ARGS[1])
n_particles = [1024, 4096, 16384, 65536, 262144, 1048576, 67108864]

for n in n_particles
    println("Pushing $n particles.")
    v_gpu = CUDA.zeros(Float64, n, 6)
    bunch = Bunch(v_gpu)
    track!(bunch, bitsring)
end
