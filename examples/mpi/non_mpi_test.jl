using Pkg
push!(LOAD_PATH, "~/repos/BeamTracking.jl")
Pkg.activate(joinpath(@__DIR__, "../.."))
Pkg.instantiate()

using BeamTracking, Beamlines, MPI, BenchmarkTools, Plots, LaTeXStrings, Unitful,
 PhysicalConstants, Random

# Read in the Electron Storage Ring of the Electron-Ion Collider
include("../../test/lattices/esr.jl") # Beamline symbol is "ring"
# Currently only Linear tracking is supported, enable it for each element
foreach(t -> t.tracking_method = Linear(), ring.line)
n_particles = 1000000 #parse(Int, ARGS[1])

start_time = time()

bunch = Bunch(n_particles)
# Track the bunch through the ESR
@time begin
	track!(bunch, ring)
end


exit(0)