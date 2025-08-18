include("../single_resonance_magnus_6.jl")
using CUDA, StaticArrays, LinearAlgebra, Plots, ColorSchemes

n_particles = 128
threads = 16
using Distributions

Y_cpu = Matrix(rand(MultivariateNormal([0; 0], 0.5I), n_particles)')
scatter(Y_cpu[:, 1], Y_cpu[:, 2])
# Y_cpu = repeat([0.0, 1.0]', n_particles)
S_cpu = repeat([0.0, 0.0, 0.0, 1.0]', n_particles)

# Transfer to GPU
Y_dev = CUDA.zeros(Float64, n_particles, 2)
S_dev = CUDA.zeros(Float64, n_particles, 4)
CUDA.copyto!(Y_dev, Y_cpu)
CUDA.copyto!(S_dev, S_cpu)

beta = 1.0
epsilon = 1 / (619 * 2 * pi)
nu0 = 4.0
dt = 1 / nu0 / 10 / 2π
delta = 0.0
sigma0 = 1.647065609e-4
tf = 1 / beta * (2π)

steps = round(Int, tf / dt)
blocks = cld(n_particles, threads)

Y_all = zeros(steps, n_particles, 2)
S_all = zeros(steps, n_particles, 4)

for i in 1:steps
    @cuda threads = threads blocks = blocks single_resonance_track_magnus(Y_dev, S_dev, nu0, sigma0, beta, epsilon, 1, dt)
    CUDA.synchronize()
    Y_cpu = Array(Y_dev)
    S_cpu = Array(S_dev)
    Y_all[i, :, :] = Y_cpu
    S_all[i, :, :] = S_cpu
end

Y_all
scatter(Y_all[:, :, 1], Y_all[:, :, 2], aspect_ratio=:equal, markersize=0.2)

proj_n(v, n) = ((sum(v .* n))) * n
proj_plane(v, n) = v - proj_n(v, n)
plot()
spin_anim = @animate for i in 1:10:steps

    tip_x = Y_all[i, :, 1] .+ S_all[i, :, 2]
    tip_y = Y_all[i, :, 2] .+ S_all[i, :, 3]
    tip_z = 0.0 .+ S_all[i, :, 4]

    spin_vectors = hcat(S_all[i, :, 2], S_all[i, :, 3], S_all[i, :, 4])
    angles = acos.(spin_vectors[:, 3]) * 2 # times 2 to get wider color spread
    cmap = cgrad(:inferno, 256)  # 256-color gradient
    colors = [RGB(get(cmap, val)) for val in angles]

    # 3D quiver plot of the vector
    p1 = quiver(Y_all[i, :, 1], Y_all[i, :, 2], [0.0], quiver=(S_all[i, :, 2], S_all[i, :, 3], S_all[i, :, 4]),
        arrow=true,
        xlim=(-3.0, 3.0),
        ylim=(-3.0, 3.0),
        zlim=(-1.2, 1.2),
        xlabel="x", ylabel="y", zlabel="z",
        title="Particle spin vector",
        label="Spin vector", camera=(30, 30), legend=true, size=(1080, 1080))
    # quiver!([Y_all[i, :, 1]], [Y_all[i, :, 2]], [0.0], quiver=([0.0], [0.0], [S_all[i, :, 4]]), label="<0, 0, 1>")

    plot!(legend=:topright)

    scatter(p1, tip_x, tip_y, tip_z, color=colors, markersize=6)
end

# plot(trail_x, trail_y)

gif(spin_anim, "spin_arrow.gif", fps=24)


# using Distributions
# I = [1.0 0.0; 0.0 1.0]
# Y_sample = rand(MultivariateNormal([0; 0], 0.5I), 1000)
# scatter!(Y_sample[1,:], Y_sample[2,:])

