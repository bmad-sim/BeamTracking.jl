using Plots

n_particles = 1024
threads = 16
using Distributions
# I = [1.0 0.0; 0.0 1.0]

beta = 62.45
epsilon = 1 / (619 * 2 * pi)
nu0 = 62.5
dt = 1e-3#parse(Float64, ARGS[1])
sigma0 = 1.647065609e-4
tf = 60000 * (2π)
zeta = (nu0 - beta) / sigma0

measurements = 100

t = range(0, tf, measurements + 1)
t = t[1:end-1]

t = collect(t)

plot(t, exp.(-5.58000000153889e-9 * t), yscale=:log10)