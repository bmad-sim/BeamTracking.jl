using LinearAlgebra, StaticArrays, CUDA, Statistics

include("src/single_resonance/single_resonance_magnus_6.jl")

"""
TODO - stop updating the spin every step and instead update the rotation quaternion
     - add 6th order integration
"""

# set variables and allocate arrays

n_particles = [1024 2048 4096]
threads = 256
using Distributions
# I = [1.0 0.0; 0.0 1.0]

beta = 62.45
epsilon = 1 / (619 * 2 * pi)
nu0 = 62.5
dt = 5e-4#parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4
orbits = 1
tf = orbits * 2pi
zeta = (nu0 - beta) / sigma0

steps_per_measurement = ceil(Int, orbits / dt)

for i = 1:length(n_particles)

    Y_cpu = Matrix(rand(MultivariateNormal([0; 0], 0.5I), n_particles[i])')

    R_cpu = repeat([1.0, 0.0, 0.0, 0.0]', n_particles[i])

    Y_dev = CUDA.zeros(Float64, n_particles[i], 2)
    R_dev = CUDA.zeros(Float64, n_particles[i], 4)
    CUDA.copyto!(Y_dev, Y_cpu)
    CUDA.copyto!(R_dev, R_cpu)


    blocks = cld(n_particles[i], threads)


    @cuda threads = threads blocks = blocks single_resonance_track_magnus_6(Y_dev, R_dev, nu0, sigma0, beta, epsilon, steps_per_measurement, dt)
    CUDA.synchronize()
end