using BeamTracking
using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                    RungeKuttaTracking, Bunch, STATE_ALIVE
using StaticArrays
using BenchmarkTools
using Random

function setup_particle(pc=1e9)
    species = Species("electron")
    mc2 = massof(species)
    p_over_q_ref = pc_to_R(species, pc)

    beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
    tilde_m = 1 / beta_gamma_0
    beta_0 = beta_gamma_0 / sqrt(1 + beta_gamma_0^2)
    charge = chargeof(species)
    p0c = R_to_pc(species, p_over_q_ref)

    return species, p_over_q_ref, beta_0, tilde_m, charge, p0c, mc2
end

function setup_solenoid_benchmark()
    species, p_over_q_ref, beta_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    L = 1.0
    ds_step = 0.01
    n_steps = 100
    gx = 0.0
    gy = 0.0

    # Solenoid field
    Bz_physical = 0.01  # Tesla
    source = MultipoleField(SA[0], SA[Bz_physical], SA[0.0])

    return bunch, beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps,
           gx, gy, source
end

function reset_bunch!(bunch)
    bunch.coords.v .= 0.0
    bunch.coords.v[1, BeamTracking.PXI] = 0.01
    bunch.coords.state[1] = STATE_ALIVE
end

# Setup
bunch, beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps,
    gx, gy, source = setup_solenoid_benchmark()

println("rk4_kernel! benchmark (1 particle)")
println("=========================================")
println("L: $L, ds_step: $ds_step, n_steps: $n_steps")
println()

# Warmup
reset_bunch!(bunch)
RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                               charge, p0c, mc2, L, ds_step, n_steps,
                               gx, gy, source)

# Benchmark
reset_bunch!(bunch)
b = @benchmark begin
    RungeKuttaTracking.rk4_kernel!(1, $bunch.coords, $beta_0, $tilde_m,
                                   $charge, $p0c, $mc2, $L, $ds_step, $n_steps,
                                   $gx, $gy, $source)
end setup=(reset_bunch!($bunch)) evals=1 seconds=10

display(b)
single_median = median(b)
println("Median time: $(single_median.time) ns")
println("Memory: $(single_median.memory) bytes")
println("Allocations: $(single_median.allocs)")
println()

# Multi-particle benchmark
println("rk4_kernel! benchmark (1000 particles)")
println("=========================================")

function setup_multi_particle(n_particles)
    species, p_over_q_ref, beta_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    rng = MersenneTwister(1234)
    bunch = Bunch(randn(rng, n_particles, 6) * 0.001,
                  p_over_q_ref=p_over_q_ref, species=species)

    L = 1.0
    ds_step = 0.01
    n_steps = 100
    gx = 0.0
    gy = 0.0

    Bz_physical = 0.01
    source = MultipoleField(SA[0], SA[Bz_physical], SA[0.0])

    return bunch, beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps,
           gx, gy, source
end

function track_all_particles!(bunch, beta_0, tilde_m, charge, p0c, mc2,
                              L, ds_step, n_steps, gx, gy,
                              source)
    n = size(bunch.coords.v, 1)
    for i in 1:n
        RungeKuttaTracking.rk4_kernel!(i, bunch.coords, beta_0, tilde_m,
                                       charge, p0c, mc2, L, ds_step, n_steps,
                                       gx, gy, source)
    end
    return nothing
end

n_particles = 1000
bunch_mp, beta_0_mp, tilde_m_mp, charge_mp, p0c_mp, mc2_mp,
    L_mp, ds_step_mp, n_steps_mp, gx_mp, gy_mp,
    source_mp = setup_multi_particle(n_particles)

# Store initial state for reset
v_init = copy(bunch_mp.coords.v)
state_init = copy(bunch_mp.coords.state)

function reset_multi!(bunch, v_init, state_init)
    bunch.coords.v .= v_init
    bunch.coords.state .= state_init
end

# Warmup
track_all_particles!(bunch_mp, beta_0_mp, tilde_m_mp, charge_mp, p0c_mp, mc2_mp,
                     L_mp, ds_step_mp, n_steps_mp, gx_mp, gy_mp,
                     source_mp)

# Benchmark
b_mp = @benchmark begin
    track_all_particles!($bunch_mp, $beta_0_mp, $tilde_m_mp, $charge_mp,
                         $p0c_mp, $mc2_mp, $L_mp, $ds_step_mp, $n_steps_mp,
                         $gx_mp, $gy_mp,
                         $source_mp)
end setup=(reset_multi!($bunch_mp, $v_init, $state_init)) evals=1 seconds=10

display(b_mp)
println()

multi_median = median(b_mp)
println("Median time: $(multi_median.time) ns")
println("Memory: $(multi_median.memory) bytes")
println("Allocations: $(multi_median.allocs)")
println("Per-particle median time: $(multi_median.time / n_particles) ns")
println("Per-step median time: $(multi_median.time / (n_particles * n_steps_mp)) ns")
println("Per-particle memory: $(multi_median.memory / n_particles) bytes")
