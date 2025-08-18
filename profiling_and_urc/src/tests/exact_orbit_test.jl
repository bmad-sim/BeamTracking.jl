using LinearAlgebra, Plots

n_particles = 1
dt = 0.1
tf = 5

epsilon = 0.5
beta = 1

steps = round(Int, tf / dt)

Y = ones(steps + 1, n_particles, 2)

for i in 2:steps + 1
    Y[i, 1] = exp(-1 * epsilon * dt) * ((cos(beta * dt) * Y[i - 1, 1]) + (-1 * sin(beta * dt) * Y[i - 1, 2])) 
    Y[i, 2] = exp(-1 * epsilon * dt) * ((sin(beta * dt) * Y[i - 1, 1]) + (cos(beta * dt) * Y[i - 1, 2]))
end

anim = @animate for n in 1:10:steps
    plot(Y[1:n, 1], Y[1:n, 2], xlim = [-1.5, 1.5], ylim = [-1.5, 1.5])
end

gif(anim, "orbit_test_animation.gif", fps = 24)


# function exact_solution(t, beta, epsilon, Y0)
#     A = [(-1 * epsilon) (-1 * beta);
#           beta (-1 * epsilon)]
    
#     return exp(A * t) * Y0
# end

# direct = exact_solution(steps * dt, beta, epsilon, [1; 1])

# println("Direct error:", norm(direct - Y[end, :]))



