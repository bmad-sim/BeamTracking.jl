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

    L = 1.0
    ds_step = 0.01
    n_steps = 100
    gx = 0.0
    gy = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   mm, kn, ks, p_over_q_ref)

    # Regression test
    solution = [0.0100005  0.01  0.0  0.0  -5.00038e-5  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Solenoid" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    L = 1.0
    ds_step = 0.01
    n_steps = 100
    gx = 0.0
    gy = 0.0
    
    # Solenoid field
    Bz_physical = 0.01  # Tesla
    Bz_normalized = Bz_physical / p_over_q_ref
    mm = SVector(0)  # Solenoid (m=0)
    kn = SVector(Bz_normalized)
    ks = SVector(0.0)

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   mm, kn, ks, p_over_q_ref)

    # In uniform B-field, particle should follow circular path
    # Total transverse momentum should be conserved
    pt2 = bunch.coords.v[1, 2]^2 + bunch.coords.v[1, 4]^2
    @test isapprox(pt2, 0.01^2, rtol=1e-4)
    # Regression test
    solution = [0.010000485056009705 0.009999955057780502 1.4991110783291216e-5 2.9980699961334158e-5 -5.000375031233078e-5 0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Dipole" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    L = 1.0
    ds_step = 0.01
    n_steps = 100
    gx = 0.0
    gy = 0.0
    
    # Dipole field
    By_physical = 0.01  # Tesla
    By_normalized = By_physical / p_over_q_ref
    mm = SVector(1)  # Dipole (m=1)
    kn = SVector(By_normalized)
    ks = SVector(0.0)

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   mm, kn, ks, p_over_q_ref)

    # Regression test
    solution = [0.011499735519796054 0.012997924579999955 0.0 0.0 -6.649432859025015e-5 0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Particle loss detection" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 1.5 # Unphysical initial momentum

    L = 1.0
    ds_step = 0.1  # 10 cm step size
    n_steps = 10
    gx = 0.0
    gy = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   mm, kn, ks, p_over_q_ref)

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

    L = 1.0
    gx = 0.0
    gy = 0.0
    
    # Empty multipole vectors for drift
    mm = SVector{0, Int}()
    kn = SVector{0, Float64}()
    ks = SVector{0, Float64}()

    # Track with different step sizes
    RungeKuttaTracking.rk4_kernel!(1, bunch1.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, 0.1, 10, gx, gy,
                                   mm, kn, ks, p_over_q_ref)
    RungeKuttaTracking.rk4_kernel!(1, bunch2.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, 0.05, 20, gx, gy,
                                   mm, kn, ks, p_over_q_ref)

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
    drift_line = Beamline([drift_ele], p_over_q_ref=p_over_q_ref, species_ref=species)

    track!(bunch, drift_line)

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
    sbend_line = Beamline([sbend_ele], p_over_q_ref=p_over_q_ref, species_ref=species)

    track!(bunch, sbend_line)

    # Regression test
    solution = [0.010000150630002367 0.009995978032305387 0.0 0.0 -0.00016899908120890584 0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
  end

  @testset "RungeKutta with different step configurations" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()

    # Test with ds_step
    drift_ds = Drift(L=1.0)
    drift_ds.tracking_method = RungeKutta(ds_step=0.1)
    line_ds = Beamline([drift_ds], p_over_q_ref=p_over_q_ref, species_ref=species)
    bunch_ds = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_ds.coords.v[1, BeamTracking.PXI] = 0.01
    track!(bunch_ds, line_ds)

    # Test with n_steps
    drift_ns = Drift(L=1.0)
    drift_ns.tracking_method = RungeKutta(n_steps=10)
    line_ns = Beamline([drift_ns], p_over_q_ref=p_over_q_ref, species_ref=species)
    bunch_ns = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_ns.coords.v[1, BeamTracking.PXI] = 0.01
    track!(bunch_ns, line_ns)

    # Both should give the same results
    @test isapprox(bunch_ds.coords.v, bunch_ns.coords.v, rtol=1e-2)
  end

  @testset "RungeKutta step selection" begin
    @test BeamTracking.find_steps(RungeKutta(ds_step=0.3), 1.0) == (4, 0.25)
    @test BeamTracking.find_steps(RungeKutta(n_steps=4), 1.0) == (4, 0.25)
  end

  @testset "Tilted reference curvature" begin
    _, p_over_q_ref, beta_0, _, tilde_m, charge, p0c, mc2 = setup_particle()
    zero_field = ntuple(_ -> 0.0, 6)

    horizontal = RungeKuttaTracking.kick_vector(
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zero_field...,
      charge, tilde_m, beta_0, 0.1, 0.0, p0c, mc2,
    )
    vertical = RungeKuttaTracking.kick_vector(
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zero_field...,
      charge, tilde_m, beta_0, 0.0, 0.1, p0c, mc2,
    )

    @test horizontal[2] ≈ 0.1
    @test horizontal[4] ≈ 0.0
    @test vertical[2] ≈ 0.0
    @test vertical[4] ≈ 0.1
  end

  @testset "RungeKutta callbacks" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    s_values = Float64[]
    ds_values = Float64[]
    function save_position!(i, coords, cur_s, cur_t_ref, cur_beta_gamma_ref,
                            last_ds_step, last_g, transforms_out!, transforms_in!)
      push!(s_values, cur_s)
      push!(ds_values, last_ds_step)
    end

    ele = Drift(L=1.0, tracking_method=RungeKutta(n_steps=4))
    line = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=species)
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species,
                  callbacks=(save_position!,))

    track!(bunch, line)

    # Three internal callbacks follow completed non-final RK steps. The common
    # tracking path supplies the fourth callback after element-exit processing.
    @test s_values ≈ [0.25, 0.5, 0.75, 1.0]
    @test ds_values ≈ fill(0.25, 4)
  end

  @testset "Per-particle reference ramp" begin
    using Beamlines

    species = Species("electron")
    E_ref = TimeDependentParam(t -> 1e9 * (1 + 1e6 * t), false)
    ele = Drift(L=1.0, tracking_method=RungeKutta(n_steps=4))
    line = Beamline([ele], E_ref=E_ref, species_ref=species)
    p_over_q_ref = line.p_over_q_ref(0.0)
    beta_gamma_ref = BeamTracking.R_to_beta_gamma(species, p_over_q_ref)
    expected_t_ref = 1.0 / BeamTracking.beta_gamma_to_v(beta_gamma_ref)
    expected_p_over_q_ref = line.p_over_q_ref(expected_t_ref)
    bunch = Bunch([0.0 0.01 0.0 0.0 0.0 0.0;
                   0.0 0.01 0.0 0.0 -1.0 0.0],
                  p_over_q_ref=p_over_q_ref, species=species)

    track!(bunch, line; ramp_update_each_particle=true)

    @test bunch.t_ref ≈ expected_t_ref
    @test bunch.p_over_q_ref ≈ expected_p_over_q_ref
    @test all(isfinite, bunch.coords.v)
  end

  @testset "Unsupported bend edge angles" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)

    entrance_edge = SBend(L=1.0, g_ref=0.1, e1=0.01,
                          tracking_method=RungeKutta())
    entrance_line = Beamline([entrance_edge], p_over_q_ref=p_over_q_ref,
                             species_ref=species)
    @test_throws ErrorException track!(bunch, entrance_line)

    exit_edge = SBend(L=1.0, g_ref=0.1, e2=0.01,
                      tracking_method=RungeKutta())
    exit_line = Beamline([exit_edge], p_over_q_ref=p_over_q_ref,
                         species_ref=species)
    @test_throws ErrorException track!(bunch, exit_line)
  end

  @testset "Zero-length elements" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()

    # Test zero-length drift should throw an error
    drift_zero = Drift(L=0.0)
    drift_zero.tracking_method = RungeKutta()
    line_zero = Beamline([drift_zero], p_over_q_ref=p_over_q_ref, species_ref=species)
    bunch_drift = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    
    @test_throws ErrorException track!(bunch_drift, line_zero)

    # Test negative length should also throw an error
    drift_negative = Drift(L=-0.1)
    drift_negative.tracking_method = RungeKutta()
    line_negative = Beamline([drift_negative], p_over_q_ref=p_over_q_ref, species_ref=species)
    bunch_negative = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    
    @test_throws ErrorException track!(bunch_negative, line_negative)
  end

end
