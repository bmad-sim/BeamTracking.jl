
@testset "FieldTracking" begin
  # Define a simple uniform electric field in x-direction
  function uniform_field(u, t, params)
    return SVector(u[2], 1.0, u[4], 0.0, u[6], 0.0)
  end
  @testset "FieldSystem!" begin

    # Test initial conditions
    du = zeros(6)
    u = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    p = (uniform_field, nothing)
    t = 0.0

    # Call field_system!
    FieldTracking.field_system!(du, u, p, t)
  end

  # Test field_track! with uniform field
  @testset "Uniform Field Tracking" begin
    # Create a single particle
    bunch = Bunch(zeros(1, 6))
    L = 1.0
    solver = Tsit5()

    # Track the particle
    FieldTracking.field_track!(1, BunchView(bunch), L, uniform_field, nothing, solver, (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false))

    # Verify final position and momentum
    @test isapprox(bunch.v[1, 1], 0.5, rtol=1e-5)  # x = x0 + 0.5*t^2
    @test isapprox(bunch.v[1, 2], 1.0, rtol=1e-5)  # px = t
  end

  # Test field_track! with multiple particles
  @testset "Multiple Particle Tracking" begin
    # Create multiple particles
    bunch = Bunch(zeros(3, 6))
    bunch.v[2, 1] = 1.0
    bunch.v[3, 2] = 1.0
    L = 1.0
    solver = RK4()
    kc = (KernelCall(FieldTracking.field_track!, (L, uniform_field, nothing, solver, (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false))),)
    # Track all particles
    BeamTracking.runkernels!(nothing, BunchView(bunch), kc)

    # Verify final positions and momenta
    @test isapprox(bunch.v[1, 1], 0.5, rtol=1e-5)
    @test isapprox(bunch.v[2, 1], 1.5, rtol=1e-5)
    @test isapprox(bunch.v[3, 1], 1.5, rtol=1e-5)
    @test isapprox(bunch.v[1, 2], 1.0, rtol=1e-5)
    @test isapprox(bunch.v[2, 2], 1.0, rtol=1e-5)
    @test isapprox(bunch.v[3, 2], 2.0, rtol=1e-5)
  end
end