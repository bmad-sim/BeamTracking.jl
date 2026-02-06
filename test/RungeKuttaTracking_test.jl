@testset "RungeKuttaTracking" begin
  using BeamTracking
  using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                      RungeKuttaTracking, Bunch, STATE_ALIVE, STATE_LOST_PZ, E_CHARGE, C_LIGHT
  using StaticArrays

  # Helper function to setup tracking parameters
  function setup_particle(pc=1e9)  # pc in eV, default corresponds to 1 GeV
    species = Species("electron")
    mc2 = massof(species)  # eV
    p_over_q_ref = pc_to_R(species, pc)

    # Calculate tracking parameters
    beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
    tilde_m = 1 / beta_gamma_0
    gamsqr_0 = 1 + beta_gamma_0^2
    beta_0 = beta_gamma_0 / sqrt(gamsqr_0)
    charge = chargeof(species)
    p0c = R_to_pc(species, p_over_q_ref)

    return species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2
  end

  @testset "RungeKutta constructor" begin
    using BeamTracking: RungeKutta

    # Test default constructor (no arguments)
    rk_default = RungeKutta()
    @test rk_default.ds_step == 0.2
    @test rk_default.n_steps == -1

    # Test constructor with ds_step only
    rk_ds = RungeKutta(ds_step=0.1)
    @test rk_ds.ds_step == 0.1
    @test rk_ds.n_steps == -1

    # Test constructor with n_steps only
    rk_ns = RungeKutta(n_steps=50)
    @test rk_ns.ds_step == -1.0
    @test rk_ns.n_steps == 50

    # Test constructor with both ds_step and n_steps (should error)
    @test_throws ErrorException RungeKutta(ds_step=0.1, n_steps=50)

    # Test constructor with explicit nothing values (should use defaults)
    rk_nothing = RungeKutta(ds_step=nothing, n_steps=nothing)
    @test rk_nothing.ds_step == 0.2
    @test rk_nothing.n_steps == -1
  end

  @testset "Pure drift" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    # Create bunch with small transverse momentum
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, 1.0)
    ds_step = 0.01
    g_bend = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   mm, kn, ks)

    # Regression test
    solution = [0.0100005  0.01  0.0  0.0  -5.00038e-5  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Solenoid" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, 1.0)
    ds_step = 0.01
    g_bend = 0.0
    
    # Solenoid field
    Bz_physical = 0.01  # Tesla
    Bz_normalized = Bz_physical / p_over_q_ref
    mm = SVector(0)  # Solenoid (m=0)
    kn = SVector(Bz_normalized)
    ks = SVector(0.0)

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   mm, kn, ks)

    # In uniform B-field, particle should follow circular path
    # Total transverse momentum should be conserved
    pt2 = bunch.coords.v[1, 2]^2 + bunch.coords.v[1, 4]^2
    @test isapprox(pt2, 0.01^2, rtol=1e-4)
    # Regression test
    solution = [0.0100005  0.01  -4.49423e-6  -8.988e-6  -5.00038e-5  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Dipole" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, 1.0)
    ds_step = 0.01
    g_bend = 0.0
    
    # Dipole field
    By_physical = 0.01  # Tesla
    By_normalized = By_physical / p_over_q_ref
    mm = SVector(1)  # Dipole (m=1)
    kn = SVector(By_normalized)
    ks = SVector(0.0)

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   mm, kn, ks)

    # Regression test
    solution = [0.00955106  0.00910124  0.0  0.0  -4.5644e-5  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Particle loss detection" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 1.5 # Unphysical initial momentum

    s_span = (0.0, 1.0)
    ds_step = 0.1  # 10 cm step size
    g_bend = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   mm, kn, ks)

    # Particle should not track
    solution = [0.0  1.5  0.0  0.0  0.0  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_LOST_PZ
  end

  @testset "Convergence test" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch1 = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch2 = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch1.coords.v[1, BeamTracking.PXI] = 0.01
    bunch2.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, 1.0)
    g_bend = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    # Track with different step sizes
    RungeKuttaTracking.rk4_kernel!(1, bunch1.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, 0.1, g_bend,
                                   mm, kn, ks)
    RungeKuttaTracking.rk4_kernel!(1, bunch2.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, 0.05, g_bend,
                                   mm, kn, ks)

    # Results should be identical
    @test isapprox(bunch1.coords.v, bunch2.coords.v, rtol=1e-2)
  end

  @testset "Beamlines integration - Drift" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    drift_ele = Drift(L=1.0)
    drift_ele.tracking_method = RungeKutta()

    track!(bunch, drift_ele)

    # Regression test
    solution = [0.0100005  0.01  0.0  0.0  -5.00038e-5  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
  end

  @testset "Beamlines integration - SBend" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    sbend_ele = SBend(L=1.0, angle=pi/132)
    sbend_ele.tracking_method = RungeKutta()

    track!(bunch, sbend_ele)

    # Regression test
    solution = [0.02548426139361667 0.040928045820272416 0.0 0.0 -0.0006063129051164828 0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
  end

  @testset "RungeKutta with different step configurations" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()

    # Test with ds_step
    drift_ds = Drift(L=1.0)
    drift_ds.tracking_method = RungeKutta(ds_step=0.1)
    bunch_ds = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_ds.coords.v[1, BeamTracking.PXI] = 0.01
    track!(bunch_ds, drift_ds)

    # Test with n_steps
    drift_ns = Drift(L=1.0)
    drift_ns.tracking_method = RungeKutta(n_steps=10)
    bunch_ns = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_ns.coords.v[1, BeamTracking.PXI] = 0.01
    track!(bunch_ns, drift_ns)

    # Both should give the same results
    @test isapprox(bunch_ds.coords.v, bunch_ns.coords.v, rtol=1e-2)
  end

  @testset "Zero-length elements" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()

    # Test zero-length drift should throw an error
    drift_zero = Drift(L=0.0)
    drift_zero.tracking_method = RungeKutta()
    bunch_drift = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    
    @test_throws ErrorException track!(bunch_drift, drift_zero)

    # Test negative length should also throw an error
    drift_negative = Drift(L=-0.1)
    drift_negative.tracking_method = RungeKutta()
    bunch_negative = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    
    @test_throws ErrorException track!(bunch_negative, drift_negative)
  end

end
