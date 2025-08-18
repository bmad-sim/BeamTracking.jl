using Plots

include("single_resonance_magnus.jl")

# set variables and allocate arrays

n_particles = 1024
threads = 16
using Distributions
# I = [1.0 0.0; 0.0 1.0]

beta = 62.45
epsilon = 1 / (619 * 2 * pi)
nu0 = 63
dt = 0.001 #parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4
tf = 3000 * (2π)
zeta = (nu0 - beta) / sigma0
measurements = 100

Y_cpu = Matrix(rand(MultivariateNormal([0; 0], 0.5I), n_particles)')
scatter(Y_cpu[:, 1], Y_cpu[:, 2])
# Y_cpu = repeat([0.0, 1.0]', n_particles)
S_cpu = repeat([0.0, 0.0, 0.0, 1.0]', n_particles)

# find ISF
function ISF(Y, zeta)
    c = 1 
    n1 = c * Y[1]
    n2 = c * Y[2]
    n3 = c * zeta
    c = sqrt(n1^2 + n2^2 + n3^2)
    n1 /= c
    n2 /= c
    n3 /= c
    return [0.0, n1, n2, n3]'
end

# set initial spin to ISF
for i in 1:n_particles
    S_cpu[i, :] = ISF(Y_cpu[i, :], zeta)
end

# Transfer to GPU
Y_dev = CUDA.zeros(Float64, n_particles, 2)
S_dev = CUDA.zeros(Float64, n_particles, 4)
CUDA.copyto!(Y_dev, Y_cpu)
CUDA.copyto!(S_dev, S_cpu)


steps = round(Int, tf / dt)
blocks = cld(n_particles, threads)
measurments = 100
steps_per_measurement = round(Int, steps / measurments)

all_Y = zeros(measurments, n_particles, 2)
all_polarizations = zeros(measurments)

for i = 1:measurments
    println("Starting measurment $i")
    if (i > 1)
    @cuda threads=threads blocks=blocks single_resonance_track_magnus(Y_dev, S_dev, nu0, sigma0, beta, epsilon, steps_per_measurement, dt)
    CUDA.synchronize()
    end
    Y_cpu = Array(Y_dev)
    all_Y[i, :, :] = Y_cpu
    S_cpu = Array(S_dev)
    # compute local ISF polarization
    step_ISF = zeros(n_particles, 4)
    polarization_sum = 0
    for i in 1:n_particles
        # step_ISF = ISF(Y_cpu[i, :], zeta)
        # polarization_sum += abs(dot(step_ISF, S_cpu[i, :]))
        polarization_sum += abs(dot([0, 0, 0, 1], S_cpu[i, :]))
    end
    all_polarizations[i] = polarization_sum / n_particles
    println(all_polarizations[i])
end

println("Polarizations: ", all_polarizations)

t = range(0, tf, measurments)
true_dist = Matrix(rand(MultivariateNormal([0; 0], 0.5I), n_particles)')
true_dist_norm = sqrt.(true_dist[:, 1].^2 + true_dist[:, 2].^2)
anim = @animate for i in 1:measurements
    x = sqrt.(all_Y[i, :, 1].^2 + all_Y[i, :, 2].^2)
    histogram(true_dist_norm, bins=100)
    histogram!(x,
        bins=100,
        c=:viridis,
        title="Measurement $i",
        xlim = (-0.0, 5.0),
        ylim = (0.0, 80.0))
    
end

gif(anim, "orbit_histogram.gif", fps=10)