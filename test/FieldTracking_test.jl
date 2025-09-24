
@testset "FieldTracking" begin
  # Define a simple uniform electric field in x-direction
  function uniform_field(u, t, params)
    return SVector(u[2], 1.0, u[4], 0.0, u[6], 0.0)
  end

  # Define a parametric field for testing field_params
  function parametric_field(u, t, params)
    E_x = params.E_x
    return SVector(u[2], E_x, u[4], 0.0, u[6], 0.0)
  end

  # Define a time-dependent field
  function time_varying_field(u, t, params)
    return SVector(u[2], t, u[4], 0.0, u[6], 0.0)
  end

  @testset "field_system!" begin
    # Test initial conditions with uniform field
    du = zeros(6)
    u = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    p = (uniform_field, nothing)
    t = 0.0

    # Call field_system!
    FieldTracking.field_system!(du, u, p, t)

    # Verify the derivatives are correctly computed
    @test du[1] ≈ 0.0  # dx/dt = px
    @test du[2] ≈ 1.0  # dpx/dt = Ex = 1.0
    @test du[3] ≈ 0.0  # dy/dt = py
    @test du[4] ≈ 0.0  # dpy/dt = Ey = 0.0
    @test du[5] ≈ 0.0  # dz/dt = pz
    @test du[6] ≈ 0.0  # dpz/dt = Ez = 0.0
  end

  @testset "field_system! with parameters" begin
    # Test field_system! with parametric field
    du = zeros(6)
    u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    params = (E_x=2.5,)
    p = (parametric_field, params)
    t = 0.0

    FieldTracking.field_system!(du, u, p, t)

    @test du[1] ≈ 0.0    # dx/dt = px
    @test du[2] ≈ 2.5    # dpx/dt = Ex = 2.5
    @test du[3] ≈ 0.0    # dy/dt = py
    @test du[4] ≈ 0.0    # dpy/dt = Ey = 0.0
    @test du[5] ≈ 0.0    # dz/dt = pz
    @test du[6] ≈ 0.0    # dpz/dt = Ez = 0.0
  end

  @testset "field_system! time-dependent" begin
    # Test field_system! with time-varying field
    du = zeros(6)
    u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    p = (time_varying_field, nothing)
    t = 1.5

    FieldTracking.field_system!(du, u, p, t)

    @test du[1] ≈ 0.0    # dx/dt = px
    @test du[2] ≈ 1.5    # dpx/dt = Ex = t = 1.5
    @test du[3] ≈ 0.0    # dy/dt = py
    @test du[4] ≈ 0.0    # dpy/dt = Ey = 0.0
    @test du[5] ≈ 0.0    # dz/dt = pz
    @test du[6] ≈ 0.0    # dpz/dt = Ez = 0.0
  end

  # Test field_track! with uniform field
  @testset "Uniform Field Tracking" begin
    # Create a single particle
    bunch = Bunch(zeros(1, 6))
    L = 1.0
    solver = Tsit5()

    # Track the particle
    FieldTracking.field_track!(1, bunch.coords, L, uniform_field, nothing, solver, (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false))

    # Verify final position and momentum
    @test isapprox(bunch.coords.v[1, 1], 0.5, rtol=1e-5)  # x = x0 + 0.5*t^2
    @test isapprox(bunch.coords.v[1, 2], 1.0, rtol=1e-5)  # px = t
  end

  # Test field_track! with parametric field
  @testset "Parametric Field Tracking" begin
    # Create a single particle
    bunch = Bunch(zeros(1, 6))
    L = 1.0
    solver = Tsit5()
    field_params = (E_x=3.0,)

    # Track the particle
    FieldTracking.field_track!(1, bunch.coords, L, parametric_field, field_params, solver, (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false))

    # Verify final position and momentum with E_x = 3.0
    @test isapprox(bunch.coords.v[1, 1], 1.5, rtol=1e-5)  # x = 0.5 * 3.0 * t^2
    @test isapprox(bunch.coords.v[1, 2], 3.0, rtol=1e-5)  # px = 3.0 * t
  end

  # Test field_track! with different solver options
  @testset "Different Solver Options" begin
    # Test with different solver parameters
    bunch1 = Bunch(zeros(1, 6))
    bunch2 = Bunch(zeros(1, 6))
    L = 1.0

    # Track with different solvers
    FieldTracking.field_track!(1, bunch1.coords, L, uniform_field, nothing, Tsit5(), (reltol=1e-6, abstol=1e-8))
    FieldTracking.field_track!(1, bunch2.coords, L, uniform_field, nothing, RK4(), (dt=0.01,))

    # Both should give similar results
    @test isapprox(bunch1.coords.v[1, 1], bunch2.coords.v[1, 1], rtol=1e-3)
    @test isapprox(bunch1.coords.v[1, 2], bunch2.coords.v[1, 2], rtol=1e-3)
  end

  # Test field_track! with initial conditions
  @testset "Different Initial Conditions" begin
    # Create particle with non-zero initial conditions
    bunch = Bunch(zeros(1, 6))
    bunch.coords.v[1, 1] = 2.0  # x0 = 2.0
    bunch.coords.v[1, 2] = 1.5  # px0 = 1.5
    bunch.coords.v[1, 3] = 0.5  # y0 = 0.5
    bunch.coords.v[1, 4] = 0.2  # py0 = 0.2

    L = 1.0
    solver = Tsit5()

    # Track the particle
    FieldTracking.field_track!(1, bunch.coords, L, uniform_field, nothing, solver, (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false))

    # Verify motion in x (with field) and y (without field)
    @test isapprox(bunch.coords.v[1, 1], 2.0 + 1.5 * 1.0 + 0.5 * 1.0^2, rtol=1e-5)  # x motion with field
    @test isapprox(bunch.coords.v[1, 2], 1.5 + 1.0, rtol=1e-5)  # px increases due to field
    @test isapprox(bunch.coords.v[1, 3], 0.5 + 0.2 * 1.0, rtol=1e-5)  # y motion without field
    @test isapprox(bunch.coords.v[1, 4], 0.2, rtol=1e-5)  # py unchanged (no field in y)
  end
end