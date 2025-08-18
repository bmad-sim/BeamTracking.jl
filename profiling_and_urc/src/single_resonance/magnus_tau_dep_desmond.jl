
using Plots, DifferentialEquations, LinearAlgebra, StaticArrays, Statistics, Distributions, SpecialFunctions, DelimitedFiles, CUDA, Peaks

include("single_resonance_magnus_6_desmond.jl")

"""
TODO - stop updating the spin every step and instead update the rotation quaternion
     - add 6th order integration
"""

# set variables and allocate arrays

n_particles = 4096
threads = 256
using Distributions
# I = [1.0 0.0; 0.0 1.0]

beta = 62.45
epsilon = 0
nu0 = parse(Float64, ARGS[1])
dt = 5e-4#parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4
orbits = 500
tf = orbits * 2pi
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
turns_per_measurement = ceil(Int, orbits / measurements)
steps_per_measurement = ceil(Int, turns_per_measurement * 2π / dt)
true_dt = turns_per_measurement * 2π / steps_per_measurement

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

true_theta = 0

for i = 1:measurements
    println(round(i / measurements * 100, digits=3), "% complete")
    polarization_sum = 0
    spin_sum = [0, 0, 0]
    if (i > 1)
        @cuda threads = threads blocks = blocks single_resonance_track_magnus_6_desmond(Y_dev, R_dev, nu0, sigma0, beta, epsilon, steps_per_measurement, true_dt, true_theta)
        # single_resonance_track_magnus_6_threads(Y_cpu, R_cpu, nu0, sigma0, beta, epsilon, steps_per_measurement, dt)
        CUDA.synchronize()
        global true_theta += true_dt * steps_per_measurement

        R_cpu = Array(R_dev)
        # println(R_cpu)
        # Y_cpu = Array(Y_dev)
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


println("Polarizations: ", all_polarizations)

t = range(0, tf, measurements + 1)
t = t[1:end-1] ./ 2pi

ENV["GKSwstype"] = "nul"
default(display_type=:inline)
plot(t[1:end], all_polarizations[1:end], xlabel="Orbit", ylabel="Polarization", label="Polarization Data")
t = collect(t)
log_polarizations = log.(all_polarizations[1:end])

using LsqFit
model1(t, param) =  param[1] * exp.(param[2] * t)
param1 = [1, -0.001]
lsq_fit = curve_fit(model1, t, all_polarizations, param1)
plot!(t, lsq_fit.param[1] * exp.(lsq_fit.param[2] * t), label="lsq_fit")

cof = pi / (619 * 2)
τ_dep_inv(ζ; cof=1) = cof * ((ζ^2 - 1) * -expintx(ζ^2) + 1)
tau_0 = 2pi^2 * 27.54^5 / 99 / 600^2 / 299792458
a = τ_dep_inv((nu0 - beta) / sigma0; cof=cof)

tau_tot_inv = tau_0 + a
tau_tot_inv_sim = tau_0 - lsq_fit.param[2]
P_eq = 92.38 * (tau_0 / tau_tot_inv)
P_eq_sim = 92.38 * (tau_0 / tau_tot_inv_sim)

println("\nDirect: ", -a)
println("Magnus solve: ", lsq_fit.param[2])
relative_error = abs(-a - lsq_fit.param[2]) / a
println("Relative error: ", relative_error)
println("Simulated equlibrium polarization: ", P_eq_sim)
plot!(t[1:end], exp.(-a * t[1:end]), label="Exact solve")

# savefig("nu_scan_desmond/polarizations_$dt" * "_$orbits" * "_$n_particles" * "_$nu0" * "_desmondISF.png")
writedlm("nu_scan_desmond/polarizations_$dt" * "_$orbits" * "_$n_particles" * "_$nu0" * "_desmondISF_0.01sine.txt", all_polarizations)

### THE FILES WITHOUT 0.01sine WERE RUN WITH 0.005 * sin((0.03 + 0.001 * randn())