@testset "FieldTracking" begin
    @testset "FieldSystem!" begin
        # Define a simple uniform electric field in x-direction
        function uniform_field(x, y, z, params)
            return SVector(1.0, 0.0, 0.0)
        end

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
        bunch = Bunch(1)
        work = zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(ele.tracking_method))
        L = 1.0
        solver = Tsit5()
        
        # Track the particle
        FieldTracking.field_track!(1, soaview(bunch), work, L, uniform_field, nothing, solver)
        
        # Verify final position and momentum
        @test isapprox(bunch.v[1,1], 0.5, rtol=1e-5)  # x = x0 + 0.5*t^2
        @test isapprox(bunch.v[1,2], 1.0, rtol=1e-5)  # px = t
    end
    
    # Test field_track! with multiple particles
    @testset "Multiple Particle Tracking" begin
        # Create multiple particles
        bunch = Bunch(zeros(3,6))
        bunch.v[2,1] = 1.0
        bunch.v[3,2] = 1.0
        work = zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(ele.tracking_method))
        L = 1.0
        solver = Tsit5()
        
        # Track all particles
        runkernel!(FieldTracking.field_track!, nothing, soaview(bunch), work, L, uniform_field, nothing, solver)
        
        # Verify final positions and momenta
        @test isapprox(bunch.v[1,1], 0.5, rtol=1e-5)
        @test isapprox(bunch.v[2,1], 1.5, rtol=1e-5)
        @test isapprox(bunch.v[3,1], 1.5, rtol=1e-5)
        @test isapprox(bunch.v[1,2], 1.0, rtol=1e-5)
        @test isapprox(bunch.v[2,2], 1.0, rtol=1e-5)
        @test isapprox(bunch.v[3,2], 2.0, rtol=1e-5)
    end
end 