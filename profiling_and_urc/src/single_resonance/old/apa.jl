include("single_resonance_magnus.jl")

n_particles = 1
threads = 1
using Distributions

# Y_cpu = Matrix(rand(MultivariateNormal([0; 0], 0.5I), 1)')
# scatter(Y_cpu[:, 1], Y_cpu[:, 2])
Y_cpu = [0.0, 1.0]
S_cpu = 0.0, 0.0, 0.0, 1.0
beta = 1.0
epsilon = 1 / (619 * 2 * pi)
nu0 = 4.0
dt = 0.01 
delta = 0.0
sigma0 = 1.647065609e-4
t0 = 0.0
tf = 2 / beta * (2π) 


Nsteps = Int(ceil(tf / dt))
Y = zeros(2, Nsteps)
Y[:,1] = Yinitial
S = zeros(4, Nsteps)
S[:, 1] = Sinitial
Y0 = Yinitial
S0 = Sinitial
for i = 2:Nsteps
    magnus_step!(Y0, S0, dt, beta, nu0, epsilon, sigma0)
    Y[:, i] = Y0
    S[:, i] = S0
    #scatter(Y[1,1:i], Y[2, 1:i],  xlims = (-1.2, 1.2), ylims = (-1.2, 1.2))
end


using Plots

Y
scatter(Y[1,:], Y[2,:], aspect_ratio = :equal, markersize = 0.2)

trail_x = Float64[]
trail_y = Float64[]
trail_z = Float64[]
spin_magnification = 1
proj_n(v, n) = ((sum(v .* n)) ) * n
proj_plane(v, n) = v - proj_n(v, n)
plot()
spin_anim = @animate for i in 1:1:Nsteps

    tip_x = Y[1, i] + S[2, i] * spin_magnification
    tip_y = Y[2, i] + S[3, i] * spin_magnification
    tip_z = 0.0     + S[4, i] * spin_magnification

    tip = proj_plane(S[2:4, i], [0, 0, 1]) 
    tip_x = tip[1] + Y[1, i] 
    tip_y = tip[2] + Y[2, i]
    # tip_z = tip[3] 

    push!(trail_x, tip_x)
    push!(trail_y, tip_y)
    push!(trail_z, tip_z)


    # 3D quiver plot of the vector
    p1 = quiver([Y[1, i]], [Y[2, i]], [0.0], quiver=([S[2, i]* spin_magnification], [S[3, i]]* spin_magnification, [S[4, i]]* spin_magnification), 
        arrow=true,
           xlim=(-3.0, 3.0),
           ylim=(-3.0, 3.0),
           zlim=(-1.2, 1.2),
           xlabel="x", ylabel="y", zlabel="z",
           title="Particle spin vector (magnified $spin_magnification times)",
           label="Spin vector", camera=(30, 30), legend=true, size=(1080, 1080))
    # quiver!([Y[1, i]], [Y[2, i]], [0.0], quiver=([0.0], [0.0], [S[4, i]]), label="<0, 0, 1>")

    plot!(legend=:topright)

    plot!(p1, trail_x, trail_y, trail_z, color=:blue)
end

plot(trail_x, trail_y)

gif(spin_anim, "spin_arrow.gif", fps=24)

# using Distributions
# I = [1.0 0.0; 0.0 1.0]
# Y_sample = rand(MultivariateNormal([0; 0], 0.5I), 1000)
# scatter!(Y_sample[1,:], Y_sample[2,:])

