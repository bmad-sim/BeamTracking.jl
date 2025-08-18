using Plots, LinearAlgebra, LaTeXStrings

# A2000 GPU
# One orbit, time step 5e-4
#12 GB memory, clock speed 562.5 Mhz
# 1024 particles took 174.77 ms on average
# 2048 particles took 174.82 ms on average
# 4096 particles took 183.81 ms on average

a2000_performance = [3930673953.02, 7863367416.72, 14953177614.08]
a2000_intensity = [1303.95, 1822.23, 2130.26]

a2000_min_performance = 2667412587.41
a2000_min_intensity = 0.01
a2000_max_performance = 56272727272.73
a2000_max_intensity = 0.21

a2000_nparticles = [1024.0, 2048, 4096]

p1 = scatter(a2000_intensity, a2000_performance, xscale=:log10, yscale=:log10,
    xlim=(0.01, 1e6), ylim=(1e6, 1e11), ylabel="Performance (FLOP/s)",
    xlabel="Arithmetic Intensity (FLOP/byte)", legend=false, markersize=10,
    zcolor=a2000_nparticles, colorbar_title="\nParticles", colorbar=true, clim=(1000, 5000), margin=5Plots.mm)
# annotate!(a2000_intensity[1] + 1e4, a2000_performance[1] - 1e9, "1024")
# annotate!(a2000_intensity[2] + 1e4, a2000_performance[2] - 1e9, "2048")
# annotate!(a2000_intensity[3] + 1e4, a2000_performance[3] - 1e9, "4096")
plot!([a2000_min_intensity, a2000_max_intensity], [a2000_min_performance, a2000_max_performance], linewidth=4, color="#000F7E")
hline!([a2000_max_performance], linewidth=4, color="#000F7E")
scatter!([a2000_min_intensity], [a2000_min_performance], markersize=5, marker=:utriangle, color="#000F7E")
scatter!([a2000_max_intensity], [a2000_max_performance], markersize=5, marker=:utriangle, color="#000F7E")
# A100-SXM4
# 10 orbits, time step 5e-4
#40gb memory, clock speed 1.09 Ghz
# 1024 2048 4096 8192 32768 65536 262144

savefig(p1, "figures/a2000_roofline.pdf")

a100_performance = [11396646639.11, 19714489274.88, 45455981501.07, 88422985762.58, 123860587056.50, 151213344550.40, 168581602217.29]
a100_intensity = [31074.11, 30656.78, 29357.33, 38789.18, 21267.12, 46962.47, 62095.75]

a100_long_performance = [9024837805.96, 74891468395.59]
a100_long_intensity = [39612.55, 43375.10]

a100_min_performance = 151466666666.67
a100_min_intensity = 0.1
a100_max_performance = 7358291666666.67
a100_max_intensity = 4.86

a100_nparticles = [1024.0, 2048, 4096, 8192, 32768, 65536, 262144]

gr()

p2 = scatter(a100_intensity, a100_performance, xscale=:log10, yscale=:log10,
    xlim=(0.1, 1e7), ylim=(1e6, 1e13), ylabel="Performance (FLOP/s)",
    xlabel="Arithmetic Intensity (FLOP/byte)",
    colorbar=true, zcolor=log2.(a100_nparticles), clim=(9, 19),
    legend=false, label="", markersize=8, colorbar_title=L"$\log_{2}$ particles", margin=5Plots.mm)

plot!([a100_min_intensity, a100_max_intensity], [a100_min_performance, a100_max_performance], linewidth=4, color="#000F7E")
hline!([a100_max_performance], linewidth=4, color="#000F7E")
scatter!([a100_min_intensity], [a100_min_performance], markersize=5, marker=:utriangle, color="#000F7E")
scatter!([a100_max_intensity], [a100_max_performance], markersize=5, marker=:utriangle, color="#000F7E")

savefig(p2, "figures/a100_roofline.pdf")