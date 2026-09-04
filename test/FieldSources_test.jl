function test_uniform_field(x, y, z, s, parameters)
  carrier = zero(x)
  return EMField(
    carrier + parameters.Ex,
    carrier + parameters.Ey,
    carrier + parameters.Ez,
    carrier + parameters.Bx,
    carrier + parameters.By,
    carrier + parameters.Bz,
  )
end

function test_parameter_free_field(x, y, z, s)
  carrier = zero(x)
  return EMField(carrier, carrier, carrier, carrier, carrier, carrier + 1)
end

@testset "Field sources" begin
  @testset "EMField" begin
    field = EMField(SA[1.0, 2.0, 3.0], SA[4.0, 5.0, 6.0])
    @test field.E == SA[1.0, 2.0, 3.0]
    @test field.B == SA[4.0, 5.0, 6.0]
    @test EMField(1, 2, 3, 4, 5, 6) ==
          EMField(SA[1, 2, 3], SA[4, 5, 6])
    @test field + field == EMField(SA[2.0, 4.0, 6.0], SA[8.0, 10.0, 12.0])
  end

  @testset "ZeroField" begin
    source = ZeroField()
    field = @inferred source(1.0, 2.0, 3.0, 4.0)
    @test field == EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 0.0])

    dual = ForwardDiff.Dual(1.0, 1.0)
    dual_field = @inferred source(dual, dual, dual, 0.0)
    @test eltype(dual_field.E) === typeof(dual)
    @test eltype(dual_field.B) === typeof(dual)
  end

  @testset "MultipoleField" begin
    solenoid = MultipoleField(SA[0], SA[1.5], SA[0.0])
    dipole = MultipoleField(SA[1], SA[2.0], SA[3.0])
    quadrupole = MultipoleField(SA[2], SA[4.0], SA[5.0])

    @test @inferred(solenoid(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.5])
    @test @inferred(dipole(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[3.0, 2.0, 0.0])
    @test @inferred(quadrupole(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[2.2, -0.7, 0.0])
  end

  @testset "FunctionalField" begin
    parameters = (Ex=1.0, Ey=2.0, Ez=3.0, Bx=4.0, By=5.0, Bz=6.0)
    source = FunctionalField(test_uniform_field, parameters)
    @test @inferred(source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[1.0, 2.0, 3.0], SA[4.0, 5.0, 6.0])

    parameter_free = FunctionalField(test_parameter_free_field)
    @test @inferred(parameter_free(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])

    time_source = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=2.0 * Time(), Bz=0.0),
    )
    lowered_time_source = BeamTracking.time_lower(time_source)
    @test BeamTracking.static_timecheck(lowered_time_source)
    evaluated_time_source = @inferred BeamTracking.teval(lowered_time_source, 0.25)
    @test @inferred(evaluated_time_source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.5, 0.0])

    batch_source = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=BatchParam([1.0, 2.0]), Bz=0.0),
    )
    lowered_batch_source = BeamTracking.batch_lower(batch_source)
    @test BeamTracking.static_batchcheck(lowered_batch_source)
    first_batch_source = @inferred BeamTracking.beval(lowered_batch_source, 1)
    second_batch_source = @inferred BeamTracking.beval(lowered_batch_source, 2)
    @test @inferred(first_batch_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 1.0, 0.0]
    @test @inferred(second_batch_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
  end

  @testset "SumField" begin
    dipole = MultipoleField(SA[1], SA[2.0], SA[0.0])
    external = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=1.0, By=0.0, Bz=3.0),
    )
    source = SumField(dipole, external)

    @test @inferred(source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[1.0, 2.0, 3.0])
    @test SumField(()) isa ZeroField
    @test SumField((ZeroField(), dipole)) === dipole

    nested = SumField(ZeroField(), SumField(dipole, external), ZeroField())
    @test nested isa SumField
    @test nested.sources == (dipole, external)

    dynamic = SumField(
      dipole,
      FunctionalField(
        test_uniform_field,
        (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=3.0 * Time(), Bz=0.0),
      ),
    )
    evaluated_dynamic = @inferred BeamTracking.teval(BeamTracking.time_lower(dynamic), 0.5)
    @test @inferred(evaluated_dynamic(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 3.5, 0.0]
  end

  @testset "MultipoleField parameters" begin
    batch_source = MultipoleField(
      SA[1],
      SA[BatchParam([2.0, 3.0])],
      SA[BatchParam(0.0)],
    )
    lowered_batch_source = BeamTracking.batch_lower(batch_source)
    @test BeamTracking.static_batchcheck(lowered_batch_source)
    first_batch_source = @inferred BeamTracking.beval(lowered_batch_source, 1)
    second_batch_source = @inferred BeamTracking.beval(lowered_batch_source, 2)
    @test @inferred(first_batch_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
    @test @inferred(second_batch_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 3.0, 0.0]

    time_source = MultipoleField(SA[1], SA[4.0 * Time()], SA[TimeDependentParam(0.0)])
    evaluated_time_source =
      @inferred BeamTracking.teval(BeamTracking.time_lower(time_source), 0.5)
    @test @inferred(evaluated_time_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
  end

  @testset "No scalar allocations" begin
    source = SumField(
      MultipoleField(SA[2], SA[4.0], SA[5.0]),
      FunctionalField(
        test_uniform_field,
        (Ex=0.0, Ey=0.0, Ez=0.0, Bx=1.0, By=0.0, Bz=3.0),
      ),
    )
    @test @ballocated($source(0.2, 0.3, 0.0, 0.0)) == 0
  end
end
