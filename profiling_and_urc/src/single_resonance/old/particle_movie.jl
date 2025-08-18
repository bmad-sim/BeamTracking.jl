using Plots, DifferentialEquations, LinearAlgebra, StaticArrays, Statistics, Distributions, SpecialFunctions, DelimitedFiles, CUDA, Peaks

include("single_resonance_magnus_6.jl")

"""
TODO - stop updating the spin every step and instead update the rotation quaternion
     - add 6th order integration
"""

# set variables and allocate arrays

n_particles = 6
threads = 6
using Distributions
# I = [1.0 0.0; 0.0 1.0]

beta = 62.45
epsilon = 1 / (619 * 2 * pi)
nu0 = 62.5#parse(Float64, ARGS[1])
dt = 1e-6#parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4
tf = 1 / (beta * 2pi)
zeta = (nu0 - beta) / sigma0

Y_cpu = Matrix(rand(MultivariateNormal([0; 0], 0.5I), n_particles)')
scatter(Y_cpu[:, 1], Y_cpu[:, 2])

R_cpu = repeat([1.0, 0.0, 0.0, 0.0]', n_particles)

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
    return [n1, n2, n3]'
end

# Transfer to GPU
Y_dev = CUDA.zeros(Float64, n_particles, 2)
R_dev = CUDA.zeros(Float64, n_particles, 4)
CUDA.copyto!(Y_dev, Y_cpu)
CUDA.copyto!(R_dev, R_cpu)

steps = round(Int, tf / dt)
blocks = cld(n_particles, threads)
measurements = 500
steps_per_measurement = round(Int, steps / measurements)


# steps_per_measurement = round(Int, steps / measurements)

spins = zeros(measurements, n_particles, 3)

# set initial spin to ISF
for i in 1:n_particles
    spins[1, i, :] = ISF(Y_cpu[i, :], zeta)
end

# set initial spin up
# for i in 1:n_particles
#     spins[1, i, :] = [0, 0, 1]
# end

# for i in 1:n_particles
#     spins[1, i, :] = randn(3)
#     spins[1, i, :] = spins[1, i, :] ./ norm(spins[1, i, :])
# end

all_polarizations = zeros(measurements)
all_orbits = zeros(measurements, n_particles, 2)
all_orbits[1, :, :] = Y_cpu

for i = 1:measurements
    println(round(i / measurements * 100, digits=3), "% complete")
    polarization_sum = 0
    spin_sum = [0, 0, 0]
    if (i > 1)
        @cuda threads = threads blocks = blocks single_resonance_track_magnus_6(Y_dev, R_dev, nu0, sigma0, beta, epsilon, steps_per_measurement, true_dt)
        # single_resonance_track_magnus_6_threads(Y_cpu, R_cpu, nu0, sigma0, beta, epsilon, steps_per_measurement, dt)
        CUDA.synchronize()

        R_cpu = Array(R_dev)
        # println(R_cpu)
        Y_cpu = Array(Y_dev)
        # S_cpu = Array(S_dev)

        # println(R_cpu)

    
        for j = 1:n_particles
            # step_ISF = ISF(Y_cpu[j, :], zeta)
            # polarization_sum += abs(dot(step_ISF, S_cpu[i, :]))
            # polarization_sum += abs(dot([0, 0, 0, 1], S_cpu[i, :]))
            # polarization_sum += norm([0, 0, 1] - [S_cpu[i, 2], S_cpu[i, 3], S_cpu[i, 4]])

            R_cpu[j, :] = R_cpu[j, :] ./ sqrt(sum(R_cpu[j, :] .^ 2))
            # r11 = 1 - 2 * R_cpu[j, 3]^2 - 2 * R_cpu[j, 4]^2
            # r12 = 2 * R_cpu[j, 2] * R_cpu[j, 3] - 2 * R_cpu[j, 1] * R_cpu[j, 4]
            # r13 = 2 * R_cpu[j, 2] * R_cpu[j, 4] - 2 * R_cpu[j, 1] * R_cpu[j, 3]
            # r21 = 2 * R_cpu[j, 2] * R_cpu[j, 3] - 2 * R_cpu[j, 1] * R_cpu[j, 4]
            # r22 = 1 - 2 * R_cpu[j, 2]^2 - 2 * R_cpu[j, 4]^2
            # r23 = 2 * R_cpu[j, 3] * R_cpu[j, 4] + 2 * R_cpu[j, 1] * R_cpu[j, 2]
            # r31 = 2 * R_cpu[j, 2] * R_cpu[j, 4] + 2 * R_cpu[j, 1] * R_cpu[j, 3]
            # r32 = 2 * R_cpu[j, 3] * R_cpu[j, 4] + 2 * R_cpu[j, 1] * R_cpu[j, 2]
            # r33 = 1 - 2 * R_cpu[j, 2]^2 - 2 * R_cpu[j, 3]^2

            # R = [r11 r12 r13; r21 r22 r23; r31 r32 r33]





            # current_spin = [0, spins[i-1, j, 1], spins[i-1, j, 2], spins[i-1, j, 3]]
            # rotate_quaternion!(R_cpu[j, :], current_spin)
            # spins[i, j, :] = current_spin[2:4]

            # spins[i-1, j, :] = spins[i-1, j, :] ./ norm(spins[i-1, j, :])

            q_r = R_cpu[j, 1]
            q_i = R_cpu[j, 2]
            q_j = R_cpu[j, 3]
            q_k = R_cpu[j, 4]

            # println([q_r, q_i, q_j, q_k])

            R = [
                1-2*(q_j^2+q_k^2) 2*(q_i*q_j-q_k*q_r) 2*(q_i*q_k+q_j*q_r);
                2*(q_i*q_j+q_k*q_r) 1-2*(q_i^2+q_k^2) 2*(q_j*q_k-q_i*q_r);
                2*(q_i*q_k-q_j*q_r) 2*(q_j*q_k+q_i*q_r) 1-2*(q_i^2+q_j^2)
            ]
            # println(eigen(R))

            spins[i, j, :] = R * spins[1, j, :]
            spin_sum += spins[i, j, :]
            all_orbits[i, j, :] = Y_cpu[j, :]
        end
    end

    # all_polarizations[i] = norm(spin_sum ./ n_particles)
    all_polarizations[i] = norm(spin_sum) / n_particles
    println(all_polarizations[i])
end

spin_sum = [0, 0, 0]
for i = 1:n_particles
    global spin_sum += spins[1, i, :]
end
all_polarizations[1] = norm(spin_sum ./ n_particles)

plot(all_orbits[:, 1, 1], all_orbits[:, 1, 2], aspect_ratio=:equal, markersize=0.2, lw=0.3)

anim = @animate for i = 1:measurements
        
    tip_x = all_orbits[i, :, 1] .+ spins[i, :, 1]
    tip_y = all_orbits[i, :, 2] .+ spins[i, :, 2]
    tip_z = 0.0 .+ spins[i, :, 3]

    spin_vectors = hcat(spins[i, :, 1], spins[i, :, 2], spins[i, :, 3])
    angles = acos.(spin_vectors[:, 3]) * 2 # times 2 to get wider color spread
    cmap = cgrad(:inferno, 256)  # 256-color gradient
    colors = [RGB(get(cmap, val)) for val in angles]

    # 3D quiver plot of the vector
    p1 = quiver(all_orbits[i, :, 1], all_orbits[i, :, 2], [0.0], quiver=(spins[i, :, 1], spins[i, :, 2], spins[i, :, 3]),
        arrow=true,
        xlim=(-3.0, 3.0),
        ylim=(-3.0, 3.0),
        zlim=(-1.2, 1.2),
        xlabel="x", ylabel="y", zlabel="z",
        title="Particle spin vector",
        label="Spin vector", camera=(30, 30), legend=true, size=(1080, 1080))
    # quiver!([all_orbits[i, :, 1]], [all_orbits[i, :, 2]], [0.0], quiver=([0.0], [0.0], [spins[i, :, 4]]), label="<0, 0, 1>")

    plot!(legend=:topright)

    scatter(p1, tip_x, tip_y, tip_z, color=colors, markersize=6)
end

gif(anim, "particle_movie.gif", fps = 10)