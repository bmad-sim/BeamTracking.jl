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
n_particles = parse(Int, ARGS[1])

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)
root = 0

if rank == 0
	start_time = time()
end

# block distribution
block_size, block_remainder = divrem(n_particles, comm_size)

if rank <= block_remainder
	block_size = block_size + 1
	offset = rank * block_size
else
	offset = block_remainder * (block_size + 1) + (rank - block_remainder) * block_size
end

rank_indices = offset+1 : offset+block_size
block_size = Int(rank_indices[end] - rank_indices[1])
println("Rank $rank has indices $rank_indices")

# collect number of particles on each rank
counts = MPI.Gather(block_size, root, comm)
if rank == root
    counts = counts .* 6 # for 6 phase space elements
end

Random.seed!(rank)
bunch = Bunch(block_size)
# Track the bunch through the ESR
track!(bunch, ring)
# Also can track! individual elements
# track!(bunch, ring; outer_particle_loop=true)

"""
Test before flatten and communication - matches Tracking_examples.jl when comm_size == 1 and n_particles == 100.
"""
# plot(
# 	scatter(bunch.v[:,1], bunch.v[:,2], label ="", xlabel = L"x", ylabel = L"p_x", markersize = 1),
# 	scatter(bunch.v[:,3], bunch.v[:,4], label ="", xlabel = L"y", ylabel = L"p_y", markersize = 1),
# 	scatter(bunch.v[:,5], bunch.v[:,6], label ="", xlabel = L"z", ylabel = L"p_z",
# 	markersize = 1),
# 	layout=(1,3), size=(600,300)
# )
# savefig("mpi_example_plot_no_flatten.png")

# vectorize bunch states for communication
flattened_v = Vector{Float64}(vec(transpose(bunch.v)))

# allocate buffer
if rank == root
	result_data = zeros(sum(counts))
	recv_buffer = VBuffer(result_data, counts)
else
	recv_buffer = VBuffer(nothing)
end

# collect data from ranks
MPI.Gatherv!(flattened_v, recv_buffer, 0, comm)
end_time = time()

# exit non-root ranks
if rank != root
	exit(0)
end
MPI.Finalize()


# decompress states vector
b0v = reshape(result_data, 6, :)'


# plot dim vs. momentum
plot(
	scatter(b0v[:,1], b0v[:,2], label ="", xlabel = L"x", ylabel = L"p_x", markersize = 1),
	scatter(b0v[:,3], b0v[:,4], label ="", xlabel = L"y", ylabel = L"p_y", markersize = 1),
	scatter(b0v[:,5], b0v[:,6], label ="", xlabel = L"z", ylabel = L"p_z",
	markersize = 1),
	layout=(1,3), size=(600,300)
)

savefig("mpi_example_plot.png")

elapsed_time = end_time - start_time

println("Run time: $elapsed_time seconds")

exit(0)