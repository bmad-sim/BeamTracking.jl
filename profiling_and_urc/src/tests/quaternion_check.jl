include("single_resonance_magnus.jl")

Y = randn(2)
beta = 62.45
epsilon = 1 / (619 * 2 * pi)
nu0 = 62.452
dt = 1e-5 #parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4

w = 1
x = (exp(sigma0 * Y[1] * dt) - exp(-sigma0 * Y[1] * dt)) / 4
y = (exp(sigma0 * Y[2] * dt) - exp(-sigma0 * Y[2] * dt)) / 4
z = (exp(nu0 * dt) - exp(-nu0 * dt)) / 4

converted = [w, x, y, z]

Y_n1 = MVector{2,Float64}(0.0, 0.0)
Y_n2 = MVector{2,Float64}(0.0, 0.0)
RS = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)

# first GL node Y
step_orbit!(Y, Y_n1, dt * (0.5 - 0.28867513459481287), beta, epsilon)
# second GL node Y
step_orbit!(Y, Y_n2, dt * (0.5 + 0.28867513459481287), beta, epsilon)

# S_dynamics!(Y_n1, Y_n2, RS, dt, nu0, sigma0)
# println("Solve: ", RS)
# println("Exact: ", converted)
# println("Ratio solve:exact: ", converted ./ RS)

# x30 = [cos(pi / 12), sin(pi / 12), 0, 0]
# y30 = [cos(pi / 12), 0, sin(pi / 12), 0]
# z30 = [cos(pi / 12), 0, 0, sin(pi / 12)]

# RS = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
# exponentiate_quaternion!([pi / 12, 0.0, 0.0], RS)
# println(RS)

# RS = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
# exponentiate_quaternion!([0.0, pi / 12, 0.0], RS)
# println(RS)

# RS = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
# exponentiate_quaternion!([0.0, 0.0, pi / 12], RS)
# println(RS)

k = zeros(3)
S = [0.0, 0.0, 1.0]
k[1] = (-1 * nu0 * S[2]) + (sigma0 * Y[2] * S[3])
k[2] = (nu0 * S[1]) - (sigma0 * Y[1] * S[3])
k[3] = (-1 * sigma0 * Y[2] * S[1]) + (sigma0 * Y[1] * S[2])
RS_rk4 = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
exponentiate_quaternion!(k, RS_rk4)
RS = MVector{4, Float64}(0.0, 0.0, 0.0, 0.0)
S_dynamics!(Y_n1, Y_n2, RS, dt, nu0, sigma0)
println("Magnus: ", RS)
println("RK4: ", RS_rk4 * dt)