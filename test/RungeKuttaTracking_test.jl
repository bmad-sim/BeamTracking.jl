# Field parameters use positional tuples, as supported by Beamlines EMFieldParams.
# Uniform-field components are ordered (Ex, Ey, Ez, Bx, By, Bz).
function rk_test_uniform_field(x, y, s, t, parameters)
  carrier = zero(x)
  return (
    carrier + parameters[1],
    carrier + parameters[2],
    carrier + parameters[3],
    carrier + parameters[4],
    carrier + parameters[5],
    carrier + parameters[6],
  )
end

function rk_test_parameter_free_field(x, y, s, t, parameters)
  carrier = zero(x)
  return (carrier, carrier, carrier, carrier, carrier, carrier + 1)
end

struct RKCustomField{T}
  strength::T
end

function (field_source::RKCustomField)(x, y, s, t, parameters=nothing)
  v = zero(x)
  return (v, v, v, v, v + field_source.strength, v)
end

@testset "RungeKuttaTracking" begin
  using BeamTracking
  using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                      Bunch, STATE_ALIVE, STATE_LOST_PZ, E_CHARGE
  using StaticArrays

  # Match the batch tests and the compiler-bug guard in beval: batched
  # explicit SIMD is unsupported on Julia < 1.11 with x86_64 CPUs.
  batch_simd_supported = !(VERSION < v"1.11" && Sys.ARCH == :x86_64)

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

  @testset "Field tuples" begin
    using BeamTracking: add_fields, scale_field, normalized_field_at
    field = (1, 2.0f0, 3.0, 4, 5.0f0, 6.0)
    evaluator = (x,y,s,t,p) -> p
    args = (0.0, 0.0, 0.0, 0.0, -0.5)
    @test (@inferred add_fields(field, field)) === (2, 4.0f0, 6.0, 8, 10.0f0, 12.0)
    @test (@inferred scale_field(field, 2)) === (2, 4.0f0, 6.0, 8, 10.0f0, 12.0)
    @test (@inferred normalized_field_at(evaluator, field, Val(true), args...)) === field
    @test (@inferred normalized_field_at((evaluator,), (field,), (Val(true),), args...)) === field
    @test (@inferred normalized_field_at((), (), (), args...)) === ntuple(_ -> 0.0, 6)
    @test (@inferred BeamTracking.zero_field(0.0f0, 0, 0, 0)) === ntuple(_ -> 0.0f0, 6)

    # Rejected derivatives must retain their own types, even when different
    # field components give different precision to the equations of motion.
    for px in (0.0f0, 2.0f0)
      rhs = @inferred BeamTracking._kick_vector(0.0f0, px, 0.0f0, 0.0f0, 0.0f0, 0.0f0,
        0.0f0, (1.0f0,2.0f0,3.0f0,4.0f0,5.0f0,6.0), 0.1f0, 0.99f0, 0.0f0, 0.0f0)
      @test rhs isa Tuple{Float32,Float64,Float32,Float64,Float32,Float32}
      px == 2 && @test all(iszero, rhs)
    end

    # Preserve the recursive (right-associated) summation order.
    large = ntuple(_ -> 1e16, 6)
    small = ntuple(_ -> 1.0, 6)
    @test normalized_field_at((evaluator, evaluator, evaluator),
      (large, map(-, large), small), (Val(true), Val(true), Val(true)), args...) ==
      map(+, large, map(+, map(-, large), small))

  end

  @testset "Field unit conventions" begin
    for T in (Float32, Float64), R in (T(-3), T(2))
      args = (T(0.02), T(-0.01), zero(T), zero(T))
      em_field = (x,y,s,t,p) -> p
      physical_parameters = T.((1,2,3,4,5,6))
      expected = em_field(args..., physical_parameters)
      converted = @inferred BeamTracking.normalized_field_at(em_field, physical_parameters, Val(false), args..., inv(R))
      direct = @inferred BeamTracking.normalized_field_at(em_field, physical_parameters ./ R, Val(true), args..., inv(R))
      @test all(isapprox.(converted, expected ./ R))
      @test all(isapprox.(converted, direct))
      @test converted isa NTuple{6,T}
      @test direct isa NTuple{6,T}
      multipole = (SA[1,2], T.(SA[0.1,0.2]), T.(SA[0,0]))
      result = @inferred BeamTracking.normalized_field_at((BeamTracking.multipole_field, em_field), (multipole, physical_parameters), (Val(true), Val(false)), args..., inv(R))
      @test all(isapprox.(result, BeamTracking.multipole_field(args..., multipole) .+ converted))
      for flags in ((Val(false), Val(false)), (Val(true), Val(true)), (Val(false), Val(true)))
        result = @inferred BeamTracking.normalized_field_at((em_field, em_field),
          (physical_parameters, reverse(physical_parameters)), flags, args..., inv(R))
        first_field = flags[1] == Val(true) ? physical_parameters : physical_parameters ./ R
        second_field = flags[2] == Val(true) ? reverse(physical_parameters) : reverse(physical_parameters) ./ R
        @test all(isapprox.(result, first_field .+ second_field))
      end
    end
  end

  @testset "On-axis electric acceleration" begin
    # Analytic on-axis electric acceleration, including the beta-dependent z term.
    m, beta0, Ez, z = 2.0, 1/sqrt(5.0), 0.03, 0.2
    rhs = BeamTracking.kick_vector(0., 0., 0., 0., z, 0., 0.,
      (0., 0., Ez, 0., 0., 0.), 1., m, beta0, 0., 0., 1., m)
    @test rhs[6] ≈ Ez / beta0
    @test rhs[5] ≈ m^2 * beta0 * Ez * z
  end

  @testset "RungeKutta constructor" begin
    using BeamTracking: RungeKutta

    # Test default constructor (no arguments)
    rk_default = RungeKutta()
    @test rk_default.ds_step == 0.2
    @test rk_default.n_steps == -1
    @test !rk_default.ibs_damping_on
    @test !rk_default.ibs_fluctuations_on

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

    @test fieldnames(RungeKutta) == (:ds_step, :n_steps, :ibs_damping_on, :ibs_fluctuations_on,
                                     :fringe_at, :multipole_fringe_on)
    @test RungeKutta(0.1, -1) == rk_ds
    @test_throws MethodError RungeKutta(em_field=rk_test_uniform_field)

  end

  @testset "Element field function parameters" begin
    using Beamlines
    species, R, _, _, _, _, _, _ = setup_particle()
    parameters = (0.0, 0.0, 0.0, 0.0, 0.004, 0.0)
    group = EMFieldParams(em_field=rk_test_uniform_field,
                                em_field_params=parameters)
    whole = Quadrupole(L=0.5, Kn1=0.1, EMFieldParams=group,
                       tracking_method=RungeKutta(n_steps=5))
    individual = Quadrupole(L=0.5, Kn1=0.1, em_field=rk_test_uniform_field,
                            em_field_params=parameters,
                            tracking_method=RungeKutta(n_steps=5))
    @test whole.EMFieldParams === group
    @test whole.em_field === rk_test_uniform_field
    @test whole.EMFieldParams ≈ individual.EMFieldParams
    @test whole.em_field_params === parameters
    @test !whole.em_field_normalized
    @test isnothing(Drift().EMFieldParams)
    @test Drift().em_field(0, 0, 0, 0) === (0, 0, 0, 0, 0, 0)
    copied = Beamlines.deepcopy_no_beamline(whole)
    @test copied.EMFieldParams ≈ group
    @test copied.EMFieldParams !== group
    copied.em_field = RKCustomField(0.003) # Rebuild when callable type changes.
    @test copied.em_field isa RKCustomField
    @test whole.em_field === rk_test_uniform_field
    copied.EMFieldParams = nothing
    @test isnothing(copied.EMFieldParams)

    initial = [0.001 0.01 -0.002 0.003 0.0 0.0]
    function tracked(ele)
      line = Beamline([ele]; p_over_q_ref=R, species_ref=species)
      bunch = Bunch(copy(initial); p_over_q_ref=R, species)
      track!(bunch, line; use_KA=false, use_explicit_SIMD=false)
      return bunch.coords.v
    end
    @test tracked(whole) ≈ tracked(individual)
    empty_group = Quadrupole(L=0.5, Kn1=0.1, EMFieldParams=EMFieldParams(),
                             tracking_method=RungeKutta(n_steps=5))
    plain = Quadrupole(L=0.5, Kn1=0.1, tracking_method=RungeKutta(n_steps=5))
    @test tracked(empty_group) ≈ tracked(plain)
    disabled_group = Quadrupole(L=0.5, Kn1=0.1, em_field=nothing,
                                tracking_method=RungeKutta(n_steps=5))
    @test tracked(disabled_group) ≈ tracked(plain)

    # Beamline children inherit the group and see updates on their parent.
    line = Beamline([whole]; p_over_q_ref=R, species_ref=species)
    @test line.line[1].em_field === rk_test_uniform_field
    whole.em_field = RKCustomField(0.004)
    whole.em_field_params = nothing
    @test line.line[1].em_field isa RKCustomField
    @test tracked(whole) ≈ tracked(individual)
    parameter_free = Quadrupole(L=0.5, Kn1=0.1,
      em_field=(x,y,s,t,p) -> begin
        @assert isnothing(p)
        RKCustomField(0.004)(x,y,s,t)
      end, tracking_method=RungeKutta(n_steps=5))
    @test tracked(parameter_free) ≈ tracked(individual)

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
    
    field_source = (BeamTracking.zero_field, nothing)

    BeamTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   field_source..., Val(false))

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
    field_source = (BeamTracking.multipole_field, (SA[0], SA[Bz_physical], SA[0.0]))

    BeamTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   field_source..., Val(false))

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
    field_source = (BeamTracking.multipole_field, (SA[1], SA[By_physical], SA[0.0]))

    BeamTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   field_source..., Val(false))

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
    
    field_source = (BeamTracking.zero_field, nothing)

    BeamTracking.rk4_kernel!(1, bunch.coords, beta_0, tilde_m,
                                   charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
                                   field_source..., Val(false))

    # Particle should not track
    solution = [0.0  1.5  0.0  0.0  0.0  0.0]
    @test isapprox(bunch.coords.v, solution, rtol=1e-6)
    @test bunch.coords.state[1] == STATE_LOST_PZ
  end

  @testset "Invalid RK momenta and rejected steps" begin
    species, R, beta0, _, m, charge, pc, mc2 = setup_particle()
    # Include healthy and previously lost particles in the same SIMD lanes.
    initial = [0.0 0.0 0.0 0.0 0.0 -1.0;
               0.0 0.0 0.0 0.0 0.0 -2.0;
               0.0 0.0 0.0 0.0 0.0  0.0;
               0.0 0.1 0.0 0.0 0.0  0.0;
               0.0 1.0 0.0 0.0 0.0  0.0;
               0.0 NaN 0.0 0.0 0.0  0.0;
               0.0 0.0 0.0 0.0 0.0  Inf;
               0.0 0.0 0.0 0.0 0.0  0.0]
    for (use_KA, use_explicit_SIMD) in ((false, false), (false, true), (true, false))
      bunch = Bunch(copy(initial); species, p_over_q_ref=R)
      bunch.coords.state[8] = BeamTracking.STATE_LOST_POS_X
      call = BeamTracking.make_kernel_call(BeamTracking.rk4_kernel!,
        (beta0, m, charge, pc, mc2, 1.0, 1.0, 1, 0.0, 0.0, BeamTracking.zero_field, nothing, Val(false)))
      BeamTracking.launch!(bunch.coords, call; use_KA, use_explicit_SIMD)
      @test bunch.coords.state == [STATE_LOST_PZ, STATE_LOST_PZ, STATE_ALIVE,
        STATE_ALIVE, STATE_LOST_PZ, STATE_LOST_PZ, STATE_LOST_PZ, BeamTracking.STATE_LOST_POS_X]
      @test isequal(bunch.coords.v[[1, 2, 5, 6, 7, 8], :], initial[[1, 2, 5, 6, 7, 8], :])
      @test bunch.coords.v[4, 1] ≈ 0.1 / sqrt(1 - 0.1^2)
    end

    # A bad intermediate stage must reject the step even if zeroing its
    # derivative would make the weighted final momentum appear valid.
    intermediate_loss = (BeamTracking.multipole_field, (SA[1], SA[1.1 * R], SA[0.0]))
    # All stage input momenta are valid here, but k4 makes the final px = -2.
    final_loss = ((x, y, s, t, R) -> begin
      v = zero(x)
      by = v + ifelse(s == 1, 12 * R, zero(R))
      (v, v, v, v, by, v)
    end, R)
    for field_source in (intermediate_loss, final_loss)
      for (use_KA, use_explicit_SIMD) in ((false, false), (false, true), (true, false))
        bunch = Bunch(zeros(8, 6); species, p_over_q_ref=R)
        call = BeamTracking.make_kernel_call(BeamTracking.rk4_kernel!,
          (beta0, m, charge, pc, mc2, 1.0, 1.0, 1, 0.0, 0.0, field_source..., Val(false)))
        BeamTracking.launch!(bunch.coords, call; use_KA, use_explicit_SIMD)
        @test all(==(STATE_LOST_PZ), bunch.coords.state)
        @test iszero(bunch.coords.v)
      end
      # Direct users of rk4_step! need the same loss handling as the kernel.
      bunch = Bunch(zeros(1, 6); species, p_over_q_ref=R)
      BeamTracking.rk4_step!(bunch.coords, 1, 0.0, 1.0, field_source..., Val(false),
        charge, m, beta0, 0.0, 0.0, pc, mc2)
      @test bunch.coords.state[1] == STATE_LOST_PZ
      @test iszero(bunch.coords.v)
    end

    # Electric deceleration may also make total momentum nonpositive midstep.
    field_source = ((x, y, s, t, p) -> begin
      v = zero(x)
      (v, v, v + 4 * pc, v, v, v)
    end, nothing)
    bunch = Bunch(zeros(1, 6); species, p_over_q_ref=R)
    BeamTracking.rk4_step!(bunch.coords, 1, 0.0, 1.0, field_source..., Val(false),
      charge, m, beta0, 0.0, 0.0, pc, mc2)
    @test bunch.coords.state[1] == STATE_LOST_PZ
    @test iszero(bunch.coords.v)
  end

  @testset "Finite field-query time for rejected momentum" begin
    species, R, beta0, _, m, charge, pc, mc2 = setup_particle()
    # The field is evaluated before a bad lane is rejected. It must not receive
    # NaN time just because another particle's longitudinal momentum is invalid.
    field_source = ((x, y, s, t, p) -> begin
      @assert all(isfinite(t))
      v = zero(x)
      (v, v, v, v, v, v + t)
    end, nothing)
    for T in (Float32, Float64), (use_KA, use_explicit_SIMD) in
        ((false, false), (false, true), (true, false))
      initial = zeros(T, 8, 6)
      initial[:, 5] .= T(0.2)
      initial[:, 6] .= T.((-1, -2, Inf, -Inf, NaN, 0, 0, 0))
      initial[6, 2] = T(0.01)
      bunch = Bunch(copy(initial); species, p_over_q_ref=R)
      bunch.coords.state[8] = BeamTracking.STATE_LOST_POS_X
      call = BeamTracking.make_kernel_call(BeamTracking.rk4_kernel!,
        (T(beta0), T(m), T(charge), T(pc), T(mc2), T(0.1), T(0.01), 10,
         zero(T), zero(T), field_source..., Val(false)))
      BeamTracking.launch!(bunch.coords, call; use_KA, use_explicit_SIMD)
      @test bunch.coords.state == [STATE_LOST_PZ, STATE_LOST_PZ, STATE_LOST_PZ,
        STATE_LOST_PZ, STATE_LOST_PZ, STATE_ALIVE, STATE_ALIVE, BeamTracking.STATE_LOST_POS_X]
      @test isequal(bunch.coords.v[[1, 2, 3, 4, 5, 8], :], initial[[1, 2, 3, 4, 5, 8], :])
      @test all(isfinite, bunch.coords.v[6:7, :])
    end
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

  @testset "Beamlines physical and normalized multipoles" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    initial = [0.001 0.01 -0.002 0.003 0.0 0.0]
    normalized_strength = 0.2
    physical_strength = normalized_strength * p_over_q_ref

    normalized_element = Quadrupole(
      L=0.5,
      Kn1=normalized_strength,
      tracking_method=RungeKutta(n_steps=5),
    )
    physical_element = Quadrupole(
      L=0.5,
      Bn1=physical_strength,
      tracking_method=RungeKutta(n_steps=5),
    )
    normalized_line = Beamline(
      [normalized_element],
      p_over_q_ref=p_over_q_ref,
      species_ref=species,
    )
    physical_line = Beamline(
      [physical_element],
      p_over_q_ref=p_over_q_ref,
      species_ref=species,
    )
    normalized_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
    physical_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)

    track!(normalized_bunch, normalized_line)
    track!(physical_bunch, physical_line)

    @test normalized_bunch.coords.v ≈ physical_bunch.coords.v
  end

  @testset "Beamlines additive field functions" begin
    species, R, _, _, _, _, _, _ = setup_particle()
    initial = [0.001 0.01 -0.002 0.003 0.0 0.0]
    external = (0.0, 0.0, 0.0, 0.0, 0.004, 0.0)
    # No fringe since the reference map is a custom field, which gets no fringe kick.
    element = Quadrupole(L=0.5, Kn1=0.2,
      em_field=rk_test_uniform_field, em_field_params=external,
      tracking_method=RungeKutta(n_steps=5, fringe_at=Fringe.NoEnd))
    # Independent complete map, without BMultipoleParams, gives the same total field.
    complete_map = (x,y,s,t,p) -> (zero(x), zero(x), zero(x),
                                        p[1]*y, p[1]*x + p[2], zero(x))
    line = Beamline([element]; p_over_q_ref=R, species_ref=species)
    function tracked(line)
      bunch = Bunch(copy(initial); p_over_q_ref=R, species)
      track!(bunch, line)
      return bunch.coords.v
    end
    previous = nothing
    for strength in (0.2, 1.0)
      element.Kn1 = strength
      reference = Drift(L=0.5, em_field=complete_map,
        em_field_params=(strength*R, external[5]),
        tracking_method=RungeKutta(n_steps=5))
      reference_line = Beamline([reference]; p_over_q_ref=R, species_ref=species)
      result = tracked(line)
      @test result ≈ tracked(reference_line)
      if !isnothing(previous)
        @test !(result ≈ previous)
      end
      previous = result
    end
  end

  @testset "Beamlines field function context" begin
    species, R, _, _, _, _, _, _ = setup_particle()
    initial = [0.001 0.01 -0.002 0.003 0.0 0.0]
    context = Context(dipole=0.003, external=0.004)
    parameters = (0.0, 0.0, 0.0, 0.0,
                  DefExpr{Float64}(c -> c.external), 0.0)
    context_element = Drift(L=0.5, Bn0=DefExpr{Float64}(c -> c.dipole),
      em_field=rk_test_uniform_field, em_field_params=parameters,
      tracking_method=RungeKutta(n_steps=5))
    context_line = Beamline([context_element]; context, p_over_q_ref=R, species_ref=species)
    for external in (0.004, 0.008)
      context.external = external
      fixed_element = Drift(L=0.5, Bn0=context.dipole,
        em_field=rk_test_uniform_field,
        em_field_params=(0.0, 0.0, 0.0, 0.0, external, 0.0),
        tracking_method=RungeKutta(n_steps=5))
      fixed_line = Beamline([fixed_element]; p_over_q_ref=R, species_ref=species)
      context_bunch = Bunch(copy(initial); p_over_q_ref=R, species)
      fixed_bunch = Bunch(copy(initial); p_over_q_ref=R, species)
      track!(context_bunch, context_line)
      track!(fixed_bunch, fixed_line)
      @test context_bunch.coords.v ≈ fixed_bunch.coords.v
    end
  end

  @testset "Functional field function tracking" begin
    species, p_over_q_ref, beta_0, _, tilde_m, charge, p0c, mc2 = setup_particle()
    context = Context(strength=0.004)
    initial = repeat([0.001 0.01 -0.002 0.003 0.0 0.0], 8, 1)
    initial[:, 5] .= range(0.0, -0.2, length=8)
    beta_gamma_ref = R_to_beta_gamma(species, p_over_q_ref)
    times = [BeamTracking.compute_time(initial[i, 5], initial[i, 6], 0.0, beta_gamma_ref)
             for i in axes(initial, 1)]
    cases = (
      (DefExpr{Float64}(c -> c.strength), fill(context.strength, 8), false),
      (DefExpr{BatchParam}(c -> BatchParam([c.strength, 2*c.strength])),
       repeat([context.strength, 2*context.strength], 4), false),
      (DefExpr{TimeDependentParam}(c -> c.strength + 1e6*Time()),
       context.strength .+ 1e6 .* times, false),
      (ForwardDiff.Dual(context.strength, 1.0), fill(context.strength, 8), true),
    )
    for (strength, expected_strengths, scalar_params) in cases
      functional = (
        (x, y, s, t, p) -> RKCustomField(p[1])(x, y, s, t),
        (strength,),
      )
      expected = similar(initial)
      for i in axes(initial, 1)
        line = Beamline([Drift(L=0.5, Bn0=expected_strengths[i], tracking_method=RungeKutta(n_steps=5, fringe_at=Fringe.NoEnd))],
                        p_over_q_ref=p_over_q_ref, species_ref=species)
        bunch = Bunch(copy(initial[i:i, :]), p_over_q_ref=p_over_q_ref, species=species)
        track!(bunch, line; use_KA=false, use_explicit_SIMD=false)
        expected[i, :] .= bunch.coords.v[1, :]
      end
      for (use_KA, use_explicit_SIMD) in ((false, false), (false, true), (true, false))
        if strength isa DefExpr{BatchParam} && use_explicit_SIMD && !batch_simd_supported
          continue
        end
        line = Beamline([Drift(L=0.5, em_field=functional[1],
                              em_field_params=functional[2],
                              tracking_method=RungeKutta(n_steps=5))],
                        context=context, p_over_q_ref=p_over_q_ref, species_ref=species)
        bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
        track!(bunch, line; scalar_params, use_KA, use_explicit_SIMD)
        @test bunch.coords.v ≈ expected
      end
    end

    field_source = (
      (x, y, s, t, p) -> RKCustomField(p[1])(x, y, s, t),
      (BatchParam([0.002, 0.004]),),
    )
    call = BeamTracking.make_kernel_call(BeamTracking.rk4_kernel!, (
      beta_0, tilde_m, charge, p0c, mc2, 0.5, 0.1, 5, 0.0, 0.0, field_source..., Val(false),
    ))
    bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
    for simd in (batch_simd_supported ? (false, true) : (false,))
      @test @ballocated(BeamTracking.launch!($bunch.coords, $call;
                       use_KA=false, use_explicit_SIMD=$simd)) == 0
    end
  end

  @testset "Float32 field function tracking" begin
    species, R, = setup_particle()
    # Check every RK stage, including SIMD lanes, not just storage.
    evaluator = (x, y, s, t, p) -> begin
      @assert eltype(x) === eltype(y) === eltype(s) === eltype(t) === Float32
      @assert eltype(p[1]) === Float32
      v = zero(x)
      (v, v, v, v, v + p[1], v)
    end
    for (use_KA, use_explicit_SIMD) in ((false, false), (false, true), (true, false))
      # Float32 batch gathers currently hit an upstream SIMD pointer-cast bug.
      # Cover SIMD with static parameters and batches on scalar/KA paths.
      strength = use_explicit_SIMD ? BatchParam(0.002) : BatchParam([0.002, 0.004])
      line = Beamline([Drift(L=0.5, Bn0=0.001, em_field=evaluator,
                              em_field_params=(strength,),
                              tracking_method=RungeKutta(n_steps=5))];
                      p_over_q_ref=R, species_ref=species)
      reference = Beamline([Drift(L=0.5, Bn0=strength + 0.001,
                                   tracking_method=RungeKutta(n_steps=5))];
                           p_over_q_ref=R, species_ref=species)
      bunch = Bunch(zeros(Float32, 16, 6); species, p_over_q_ref=R)
      expected = Bunch(zeros(16, 6); species, p_over_q_ref=R)
      track!(expected, reference; use_KA=false, use_explicit_SIMD=false)
      track!(bunch, line; use_KA, use_explicit_SIMD)
      # Near-zero z includes cancellation of order-one terms in Float32.
      @test all(isapprox.(bunch.coords.v, expected.coords.v; rtol=1e-5, atol=eps(Float32)))
    end
  end

  @testset "Batch field function tracking" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    batch_fields = [0.002, 0.004]
    initial_particle = [0.001 0.01 -0.002 0.003 0.0 0.0]
    initial = repeat(initial_particle, 8, 1)
    field_source = (
      rk_test_uniform_field,
      (0.0, 0.0, 0.0, 0.0, BatchParam(batch_fields), 0.0),
    )
    element = Drift(
      L=0.5,
      em_field=field_source[1], em_field_params=field_source[2],
      tracking_method=RungeKutta(n_steps=5),
    )
    line = Beamline([element], p_over_q_ref=p_over_q_ref, species_ref=species)
    simd_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
    ka_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)

    track!(simd_bunch, line; use_KA=false, use_explicit_SIMD=batch_simd_supported)
    track!(ka_bunch, line; use_KA=true, use_explicit_SIMD=false)

    expected = similar(initial)
    for i in axes(initial, 1)
      fixed_field_source = (
        rk_test_uniform_field,
        (
          0.0,
          0.0,
          0.0,
          0.0,
          batch_fields[mod1(i, length(batch_fields))],
          0.0,
        ),
      )
      fixed_element = Drift(
        L=0.5,
        em_field=fixed_field_source[1], em_field_params=fixed_field_source[2],
        tracking_method=RungeKutta(n_steps=5),
      )
      fixed_line = Beamline(
        [fixed_element],
        p_over_q_ref=p_over_q_ref,
        species_ref=species,
      )
      fixed_bunch = Bunch(
        copy(initial[i:i, :]),
        p_over_q_ref=p_over_q_ref,
        species=species,
      )
      track!(fixed_bunch, fixed_line; use_KA=false, use_explicit_SIMD=false)
      expected[i, :] .= fixed_bunch.coords.v[1, :]
    end

    @test simd_bunch.coords.v ≈ expected
    @test ka_bunch.coords.v ≈ expected

    # The same offset must reach field-function tuples through both the
    # scalar and SIMD/KA launch paths, with cyclic batch selection preserved.
    for batch_start in (2, length(batch_fields) + 2)
      shifted = expected[[mod1(i + batch_start - 1, length(batch_fields))
                          for i in axes(initial, 1)], :]
      for (use_KA, use_explicit_SIMD) in ((false, false), (false, batch_simd_supported), (true, false))
        bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
        track!(bunch, line; batch_start, use_KA, use_explicit_SIMD)
        @test bunch.coords.v ≈ shifted
      end
    end
  end

  @testset "Scalarized field function tracking" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    initial = [0.001 0.01 -0.002 0.003 0.0 0.0]
    dual_field = ForwardDiff.Dual(0.004, 1.0)
    field_source = (
      rk_test_uniform_field,
      (0.0, 0.0, 0.0, 0.0, dual_field, 0.0),
    )
    fixed_field_source = (
      rk_test_uniform_field,
      (0.0, 0.0, 0.0, 0.0, ForwardDiff.value(dual_field), 0.0),
    )
    element = Drift(
      L=0.5,
      em_field=field_source[1], em_field_params=field_source[2],
      tracking_method=RungeKutta(n_steps=5),
    )
    fixed_element = Drift(
      L=0.5,
      em_field=fixed_field_source[1], em_field_params=fixed_field_source[2],
      tracking_method=RungeKutta(n_steps=5),
    )
    line = Beamline([element], p_over_q_ref=p_over_q_ref, species_ref=species)
    fixed_line = Beamline(
      [fixed_element],
      p_over_q_ref=p_over_q_ref,
      species_ref=species,
    )
    bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)
    fixed_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)

    track!(bunch, line; scalar_params=true)
    track!(fixed_bunch, fixed_line)

    @test bunch.coords.v ≈ fixed_bunch.coords.v
  end

  @testset "Time-dependent field function tracking" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    beta_gamma_ref = R_to_beta_gamma(species, p_over_q_ref)
    initial = [
      0.001 0.01 -0.002 0.003  0.00  0.0
      0.001 0.01 -0.002 0.003 -0.05  0.0
      0.001 0.01 -0.002 0.003 -0.10  0.0
      0.001 0.01 -0.002 0.003 -0.15  0.0
    ]
    field_at_time = 0.002 + 1.0e6 * Time()
    field_source = (
      rk_test_uniform_field,
      (0.0, 0.0, 0.0, 0.0, field_at_time, 0.0),
    )
    element = Drift(
      L=0.5,
      em_field=field_source[1], em_field_params=field_source[2],
      tracking_method=RungeKutta(n_steps=5),
    )
    line = Beamline([element], p_over_q_ref=p_over_q_ref, species_ref=species)
    dynamic_bunch = Bunch(copy(initial), p_over_q_ref=p_over_q_ref, species=species)

    track!(dynamic_bunch, line; use_KA=false, use_explicit_SIMD=true)

    expected = similar(initial)
    for i in axes(initial, 1)
      particle_time = BeamTracking.compute_time(
        initial[i, BeamTracking.ZI],
        initial[i, BeamTracking.PZI],
        0.0,
        beta_gamma_ref,
      )
      fixed_field_source = (
        rk_test_uniform_field,
        (
          0.0,
          0.0,
          0.0,
          0.0,
          field_at_time(particle_time),
          0.0,
        ),
      )
      fixed_element = Drift(
        L=0.5,
        em_field=fixed_field_source[1], em_field_params=fixed_field_source[2],
        tracking_method=RungeKutta(n_steps=5),
      )
      fixed_line = Beamline(
        [fixed_element],
        p_over_q_ref=p_over_q_ref,
        species_ref=species,
      )
      fixed_bunch = Bunch(
        copy(initial[i:i, :]),
        p_over_q_ref=p_over_q_ref,
        species=species,
      )
      track!(fixed_bunch, fixed_line; use_KA=false, use_explicit_SIMD=false)
      expected[i, :] .= fixed_bunch.coords.v[1, :]
    end

    @test dynamic_bunch.coords.v ≈ expected
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

    horizontal = BeamTracking.kick_vector(
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zero_field...,
      charge, tilde_m, beta_0, 0.1, 0.0, p0c, mc2,
    )
    vertical = BeamTracking.kick_vector(
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

  @testset "Bend edge angles without the fringe at that end" begin
    using Beamlines

    species, p_over_q_ref, _, _, _, _, _, _ = setup_particle()
    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)

    entrance_edge = SBend(L=1.0, g_ref=0.1, e1=0.01,
                          tracking_method=RungeKutta(fringe_at=Fringe.ExitEnd))
    entrance_line = Beamline([entrance_edge], p_over_q_ref=p_over_q_ref,
                             species_ref=species)
    @test_throws ErrorException track!(bunch, entrance_line)

    exit_edge = SBend(L=1.0, g_ref=0.1, e2=0.01,
                      tracking_method=RungeKutta(fringe_at=Fringe.EntranceEnd))
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

@testset "RK physical and normalized functional fields" begin
  for species in (Species("electron"), Species("proton")), T in (Float32, Float64)
    R = T(chargeof(species) * 3)
    params = T.((1e4, -2e4, 3e4, 0.001, -0.002, 0.003))
    em_field = (x,y,s,t,p) -> p
    results = map(((params, false), (params ./ R, true))) do (parameters, normalized)
      ele = Drift(L=0.2, em_field=em_field, em_field_params=parameters,
        em_field_normalized=normalized, tracking_method=RungeKutta(n_steps=10))
      line = Beamline([ele]; p_over_q_ref=R, species_ref=species)
      bunch = Bunch(T.([0.001 0.002 -0.003 0.001 0.0 0.01]); p_over_q_ref=R, species=species)
      track!(bunch, line)
      bunch.coords
    end
    @test results[1].v ≈ results[2].v
    @test results[1].state == results[2].state
  end
end
