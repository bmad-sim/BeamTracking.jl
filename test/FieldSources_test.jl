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

struct FieldSourceTestAdaptor end

function BeamTracking.Adapt.adapt_storage(::FieldSourceTestAdaptor, values::Vector)
  return SVector{length(values)}(values)
end

struct FieldTestParameters{T}
  strength::T
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

    @test_throws MethodError MultipoleField(SA[1], MVector(0.01), SA[0.0])
    @test_throws MethodError MultipoleField(MVector(1), SA[0.01], SA[0.0])
    @test_throws MethodError MultipoleField(SA[1], SA[0.01], MVector(0.0))

    @test @inferred(solenoid(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.5])
    @test @inferred(dipole(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[3.0, 2.0, 0.0])
    @test @inferred(quadrupole(0.2, 0.3, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[2.2, -0.7, 0.0])

    simd_x = SIMD.Vec{2,Float64}((0.2, 0.4))
    simd_y = SIMD.Vec{2,Float64}((0.3, 0.1))
    simd_field = @inferred quadrupole(simd_x, simd_y, zero(simd_x), zero(simd_x))
    @test all(isapprox.(Tuple(simd_field.B[1]), (2.2, 2.4)))
    @test all(isapprox.(Tuple(simd_field.B[2]), (-0.7, 1.1)))

    descriptor = Descriptor(2, 1)
    tpsa_x, tpsa_y = @vars(descriptor)
    tpsa_field = @inferred quadrupole(tpsa_x, tpsa_y, zero(tpsa_x), zero(tpsa_x))
    @test GTPSA.jacobian(collect(tpsa_field.B[1:2])) ≈ [5.0 4.0; 4.0 -5.0]
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
    @test_opt source(0.2, 0.3, 0.0, 0.0)
    @test @ballocated($source(0.2, 0.3, 0.0, 0.0)) == 0
  end

  @testset "Field-source preparation" begin
    ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
    context = Beamlines.Context(strength=0.25)
    parameters = (
      strength=Beamlines.DefExpr{Float64}(c -> 2 * c.strength),
      constants=(label="map", count=2),
    )
    source = FunctionalField((x, y, z, s, p) -> begin
      v = zero(x)
      return EMField(v, v, v, v, v + p.strength, v)
    end, parameters)
    prepared = @inferred ext.unpack_field_source(source, context)
    @test prepared.evaluator === source.evaluator
    @test prepared.parameters.strength == 0.5
    @test prepared.parameters.constants === parameters.constants
    @test source.parameters.strength isa Beamlines.DefExpr

    opaque_parameters = FieldTestParameters(parameters.strength)
    opaque_source = FunctionalField(source.evaluator, opaque_parameters)
    @test ext.unpack_field_source(opaque_source, context).parameters === opaque_parameters
    @test BeamTracking.batch_lower(
      FunctionalField(source.evaluator, FieldTestParameters(BatchParam([0.5, 1.5]))),
    ).parameters isa FieldTestParameters{BatchParam}

    context.strength = 0.75
    @test ext.unpack_field_source(source, context).parameters.strength == 1.5

    dual_source = FunctionalField(source.evaluator, (
      strength=ForwardDiff.Dual(0.5, 1.0),
      constants=parameters.constants,
    ))
    scalar_source = @inferred ext.scalarize_field_source(dual_source)
    @test scalar_source.evaluator === source.evaluator
    @test scalar_source.parameters.strength === 0.5

    dynamic = FunctionalField(source.evaluator, (
      strength=BatchParam([0.5, 1.5]),
      constants=parameters.constants,
    ))
    lowered = BeamTracking.batch_lower(dynamic)
    @test @inferred(BeamTracking.static_batchcheck(lowered))
    selected = @inferred BeamTracking.beval(lowered, 2)
    @test selected(0.0, 0.0, 0.0, 0.0).B[2] == 1.5
    @test_opt BeamTracking.beval(lowered, 2)
    @test @ballocated(BeamTracking.beval($lowered, 2)) == 0

    timed = FunctionalField(source.evaluator, (
      strength=2 * Time(),
      constants=parameters.constants,
    ))
    lowered_time = BeamTracking.time_lower(timed)
    @test @inferred(BeamTracking.static_timecheck(lowered_time))
    evaluated = @inferred BeamTracking.teval(lowered_time, 0.25)
    @test evaluated(0.0, 0.0, 0.0, 0.0).B[2] == 0.5
    @test @ballocated(BeamTracking.teval($lowered_time, 0.25)) == 0

    opaque = let value = parameters.strength
      () -> value
    end
    @test ext.unpack_field_source(opaque, context) === opaque
  end

  @testset "Adaptation" begin
    source = SumField(
      MultipoleField(SA[1], SA[0.01], SA[0.0]),
      FunctionalField(test_uniform_field, (field_map=[1.0, 2.0, 3.0],)),
    )
    adapted = BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), source)

    @test adapted isa SumField
    @test adapted.sources[1] isa MultipoleField
    @test adapted.sources[2] isa FunctionalField
    @test adapted.sources[2].parameters.field_map == SA[1.0, 2.0, 3.0]
  end
end
