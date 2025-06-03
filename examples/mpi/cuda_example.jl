using Pkg
push!(LOAD_PATH, "~/repos/BeamTracking.jl")
Pkg.activate(joinpath(@__DIR__, "../.."))
Pkg.instantiate()

using BeamTracking, Beamlines, Random, CUDA

# Read in the Electron Storage Ring of the Electron-Ion Collider
include("../../test/lattices/esr.jl") # Beamline symbol is "ring"
# Currently only Linear tracking is supported, enable it for each element
foreach(t -> t.tracking_method = Linear(), ring.line)
n_particles = 10000 #parse(Int, ARGS[1])
bitsring = BitsBeamline(ring)
CUDA.@profile begin
    CUDA.@time begin
    work = CuArray{Float64}(undef, (n_particles, 1))
    v_gpu = CUDA.zeros(Float64, n_particles, 6)
    bunch = Bunch(v_gpu)
    track!(bunch, bitsring; work)
    end
end
# Track the bunch through the ESR
#@time begin
#	track!(bunch, ring)
#end


#exit(0)