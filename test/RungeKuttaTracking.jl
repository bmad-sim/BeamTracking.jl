@testset "RungeKuttaTracking" begin
  # Define a simple uniform electric field in x-direction
  function uniform_field(u, t, params)
    return SVector(u[2], 1.0, u[4], 0.0, u[6], 0.0)
  end

  # Define a time-dependent field for testing
  function time_varying_field(u, t, params)
    return SVector(u[2], t, u[4], 0.0, u[6], 0.0)
  end

  # Define a parametric field for testing
  function parametric_field(u, t, params)
    E_x = params.E_x
    return SVector(u[2], E_x, u[4], 0.0, u[6], 0.0)
  end

  @testset "rk4_step!" begin
    # Test single RK4 step with uniform field
    u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    t = 0.0
    h = 0.1
    params = nothing

    # Perform one RK4 step
    RungeKuttaTracking.rk4_step!(u, t, h, uniform_field, params)

    # Verify results (analytical solution for uniform field)
    expected_x = 0.5 * h^2  # x = 0.5 * t^2
    expected_px = h         # px = t
    
    @test isapprox(u[1], expected_x, rtol=1e-10)
    @test isapprox(u[2], expected_px, rtol=1e-10)
    @test u[3] ≈ 0.0
    @test u[4] ≈ 0.0
    @test u[5] ≈ 0.0
    @test u[6] ≈ 0.0
  end

  @testset "rk4_step! with parameters" begin
    # Test RK4 step with parametric field
    u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    t = 0.0
    h = 0.1
    params = (E_x = 2.0,)

    RungeKuttaTracking.rk4_step!(u, t, h, parametric_field, params)

    # With E_x = 2.0, acceleration is double
    expected_x = 0.5 * 2.0 * h^2  # x = 0.5 * 2.0 * t^2
    expected_px = 2.0 * h         # px = 2.0 * t
    
    @test isapprox(u[1], expected_x, rtol=1e-10)
    @test isapprox(u[2], expected_px, rtol=1e-10)
  end

  @testset "rk4_track! single particle" begin
    # Create a single particle bunch
    bunch = Bunch(zeros(1, 6))
    t_span = (0.0, 1.0)
    n_steps = 100
    params = nothing

    # Track the particle
    RungeKuttaTracking.rk4_track!(1, BunchView(bunch), t_span, uniform_field, params, n_steps)

    # Verify final position and momentum (analytical solution)
    @test isapprox(bunch.v[1, 1], 0.5, rtol=1e-5)  # x = 0.5 * t^2 at t=1
    @test isapprox(bunch.v[1, 2], 1.0, rtol=1e-5)  # px = t at t=1
    @test bunch.v[1, 3] ≈ 0.0
    @test bunch.v[1, 4] ≈ 0.0
    @test bunch.v[1, 5] ≈ 0.0
    @test bunch.v[1, 6] ≈ 0.0
  end

  @testset "rk4_track! with different initial conditions" begin
    # Create a particle with initial position and momentum
    bunch = Bunch(zeros(1, 6))
    bunch.v[1, 1] = 1.0  # x0 = 1.0
    bunch.v[1, 2] = 0.5  # px0 = 0.5
    
    t_span = (0.0, 1.0)
    n_steps = 100
    params = nothing

    RungeKuttaTracking.rk4_track!(1, BunchView(bunch), t_span, uniform_field, params, n_steps)

    # Analytical solution: x = x0 + px0*t + 0.5*E*t^2, px = px0 + E*t
    expected_x = 1.0 + 0.5 * 1.0 + 0.5 * 1.0 * 1.0^2  # 2.0
    expected_px = 0.5 + 1.0 * 1.0  # 1.5
    
    @test isapprox(bunch.v[1, 1], expected_x, rtol=1e-5)
    @test isapprox(bunch.v[1, 2], expected_px, rtol=1e-5)
  end

  @testset "rk4_track! with parameters" begin
    # Test tracking with parametric field
    bunch = Bunch(zeros(1, 6))
    t_span = (0.0, 1.0)
    n_steps = 100
    params = (E_x = 3.0,)

    RungeKuttaTracking.rk4_track!(1, BunchView(bunch), t_span, parametric_field, params, n_steps)

    # With E_x = 3.0
    expected_x = 0.5 * 3.0 * 1.0^2  # 1.5
    expected_px = 3.0 * 1.0         # 3.0
    
    @test isapprox(bunch.v[1, 1], expected_x, rtol=1e-5)
    @test isapprox(bunch.v[1, 2], expected_px, rtol=1e-5)
  end

  @testset "rk4_track! multiple particles" begin
    # Create multiple particles with different initial conditions
    bunch = Bunch(zeros(3, 6))
    bunch.v[2, 1] = 1.0  # Different initial x
    bunch.v[3, 2] = 1.0  # Different initial px
    
    t_span = (0.0, 1.0)
    n_steps = 50
    params = nothing

    # Create kernel call to track all particles
    kc = (KernelCall(RungeKuttaTracking.rk4_track!, (t_span, uniform_field, params, n_steps)),)
    BeamTracking.runkernels!(nothing, BunchView(bunch), kc)

    # Verify results for each particle
    @test isapprox(bunch.v[1, 1], 0.5, rtol=1e-5)     # Particle 1: x0=0, px0=0
    @test isapprox(bunch.v[1, 2], 1.0, rtol=1e-5)
    
    @test isapprox(bunch.v[2, 1], 1.5, rtol=1e-5)     # Particle 2: x0=1, px0=0
    @test isapprox(bunch.v[2, 2], 1.0, rtol=1e-5)
    
    @test isapprox(bunch.v[3, 1], 1.5, rtol=1e-5)     # Particle 3: x0=0, px0=1
    @test isapprox(bunch.v[3, 2], 2.0, rtol=1e-5)
  end

  @testset "rk4_track! different step sizes" begin
    # Test convergence with different step sizes
    bunch1 = Bunch(zeros(1, 6))
    bunch2 = Bunch(zeros(1, 6))
    
    t_span = (0.0, 1.0)
    params = nothing

    # Track with different step sizes
    RungeKuttaTracking.rk4_track!(1, BunchView(bunch1), t_span, uniform_field, params, 50)
    RungeKuttaTracking.rk4_track!(1, BunchView(bunch2), t_span, uniform_field, params, 200)

    # Results should be similar but more accurate with smaller steps
    @test isapprox(bunch1.v[1, 1], bunch2.v[1, 1], rtol=1e-3)
    @test isapprox(bunch1.v[1, 2], bunch2.v[1, 2], rtol=1e-3)
    
    # Both should be close to analytical solution
    @test isapprox(bunch2.v[1, 1], 0.5, rtol=1e-6)
    @test isapprox(bunch2.v[1, 2], 1.0, rtol=1e-6)
  end

  @testset "rk4_track! time-varying field" begin
    # Test with time-dependent field
    bunch = Bunch(zeros(1, 6))
    t_span = (0.0, 2.0)
    n_steps = 200
    params = nothing

    RungeKuttaTracking.rk4_track!(1, BunchView(bunch), t_span, time_varying_field, params, n_steps)

    # With time-varying field E(t) = t, analytical solution is more complex
    # Approximate verification that particle moved and gained momentum
    @test bunch.v[1, 1] > 1.0  # Should have moved in x
    @test bunch.v[1, 2] > 1.0  # Should have gained momentum
  end
end