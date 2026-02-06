using BeamTracking
using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                    RungeKuttaTracking, Bunch, STATE_ALIVE
using StaticArrays
using BenchmarkTools

function setup_particle(pc=1e9)
    species = Species("electron")
    mc2 = massof(species)
    R_ref = pc_to_R(species, pc)

    beta_gamma_0 = R_to_beta_gamma(species, R_ref)
    tilde_m = 1 / beta_gamma_0
    gamsqr_0 = 1 + beta_gamma_0^2
    beta_0 = beta_gamma_0 / sqrt(gamsqr_0)
    charge = chargeof(species)
    p0c = R_to_pc(species, R_ref)

    return species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2
end

function setup_solenoid_benchmark()
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, 1.0)
    ds_step = 0.01
    g_bend = 0.0

    # Solenoid field
    Bz_physical = 0.01  # Tesla
    Bz_normalized = Bz_physical / R_ref
    mm = SVector(0)
    kn = SVector(Bz_normalized)
    ks = SVector(0.0)

    return bunch, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2, s_span, ds_step, g_bend, mm, kn, ks
end

function reset_bunch!(bunch)
    bunch.coords.v .= 0.0
    bunch.coords.v[1, BeamTracking.PXI] = 0.01
    bunch.coords.state[1] = STATE_ALIVE
end

# Setup
bunch, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2, s_span, ds_step, g_bend, mm, kn, ks = setup_solenoid_benchmark()

println("rk4_kernel! benchmark (1 particle)")
println("=========================================")
println("s_span: $s_span, ds_step: $ds_step")
println("n_steps: $(Int(ceil((s_span[2] - s_span[1]) / ds_step)))")
println()

# Warmup
reset_bunch!(bunch)
RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                               charge, p0c, mc2, s_span, ds_step, g_bend,
                               mm, kn, ks)

# Benchmark
reset_bunch!(bunch)
b = @benchmark begin
    RungeKuttaTracking.rk4_kernel!(1, $bunch.coords, $beta_0, $gamsqr_0, $tilde_m,
                                   $charge, $p0c, $mc2, $s_span, $ds_step, $g_bend,
                                   $mm, $kn, $ks)
end setup=(reset_bunch!($bunch))

display(b)
println()

# Multi-particle benchmark
println("rk4_kernel! benchmark (1000 particles)")
println("=========================================")

function setup_multi_particle(n_particles)
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(randn(n_particles, 6) * 0.001, R_ref=R_ref, species=species)

    s_span = (0.0, 1.0)
    ds_step = 0.01
    g_bend = 0.0

    Bz_physical = 0.01
    Bz_normalized = Bz_physical / R_ref
    mm = SVector(0)
    kn = SVector(Bz_normalized)
    ks = SVector(0.0)

    return bunch, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2, s_span, ds_step, g_bend, mm, kn, ks
end

function track_all_particles!(bunch, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2,
                              s_span, ds_step, g_bend, mm, kn, ks)
    n = size(bunch.coords.v, 1)
    for i in 1:n
        RungeKuttaTracking.rk4_kernel!(i, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                       charge, p0c, mc2, s_span, ds_step, g_bend,
                                       mm, kn, ks)
    end
end

n_particles = 1000
bunch_mp, beta_0_mp, gamsqr_0_mp, tilde_m_mp, charge_mp, p0c_mp, mc2_mp,
    s_span_mp, ds_step_mp, g_bend_mp, mm_mp, kn_mp, ks_mp = setup_multi_particle(n_particles)

# Store initial state for reset
v_init = copy(bunch_mp.coords.v)
state_init = copy(bunch_mp.coords.state)

function reset_multi!(bunch, v_init, state_init)
    bunch.coords.v .= v_init
    bunch.coords.state .= state_init
end

# Warmup
track_all_particles!(bunch_mp, beta_0_mp, gamsqr_0_mp, tilde_m_mp, charge_mp, p0c_mp, mc2_mp,
                     s_span_mp, ds_step_mp, g_bend_mp, mm_mp, kn_mp, ks_mp)

# Benchmark
b_mp = @benchmark begin
    track_all_particles!($bunch_mp, $beta_0_mp, $gamsqr_0_mp, $tilde_m_mp, $charge_mp,
                         $p0c_mp, $mc2_mp, $s_span_mp, $ds_step_mp, $g_bend_mp,
                         $mm_mp, $kn_mp, $ks_mp)
end setup=(reset_multi!($bunch_mp, $v_init, $state_init))

display(b_mp)
println()

# Per-particle timing
median_time_ns = median(b_mp).time
println("\nPer-particle median time: $(median_time_ns / n_particles) ns")

reset_multi!(bunch_mp, v_init, state_init)
num_allocs_mp = @allocated track_all_particles!(bunch_mp, beta_0_mp, gamsqr_0_mp, tilde_m_mp, charge_mp,
                                                 p0c_mp, mc2_mp, s_span_mp, ds_step_mp, g_bend_mp,
                                                 mm_mp, kn_mp, ks_mp)
println("Total allocations for $n_particles particles: $num_allocs_mp bytes")
println("Per-particle allocations: $(num_allocs_mp / n_particles) bytes")
