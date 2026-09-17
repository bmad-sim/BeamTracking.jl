function test_uniform_field(x, y, s, t, parameters)
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

function test_parameter_free_field(x, y, s, t)
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
    field_source = ZeroField()
    field = @inferred field_source(1.0, 2.0, 3.0, 4.0)
    @test field == EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 0.0])

    dual = ForwardDiff.Dual(1.0, 1.0)
    dual_field = @inferred field_source(dual, dual, dual, 0.0)
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
    field_source = FunctionalField(test_uniform_field, parameters)
    @test @inferred(field_source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[1.0, 2.0, 3.0], SA[4.0, 5.0, 6.0])

    parameter_free = FunctionalField(test_parameter_free_field)
    @test @inferred(parameter_free(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])

    time_field_source = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=2.0 * Time(), Bz=0.0),
    )
    lowered_time_field_source = BeamTracking.time_lower(time_field_source)
    @test BeamTracking.static_timecheck(lowered_time_field_source)
    evaluated_time_field_source = @inferred BeamTracking.teval(lowered_time_field_source, 0.25)
    @test @inferred(evaluated_time_field_source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.5, 0.0])

    batch_field_source = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=BatchParam([1.0, 2.0]), Bz=0.0),
    )
    lowered_batch_field_source = BeamTracking.batch_lower(batch_field_source)
    @test BeamTracking.static_batchcheck(lowered_batch_field_source)
    first_batch_field_source = @inferred BeamTracking.beval(lowered_batch_field_source, 1)
    second_batch_field_source = @inferred BeamTracking.beval(lowered_batch_field_source, 2)
    @test @inferred(first_batch_field_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 1.0, 0.0]
    @test @inferred(second_batch_field_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
  end

  @testset "SumField" begin
    dipole = MultipoleField(SA[1], SA[2.0], SA[0.0])
    external = FunctionalField(
      test_uniform_field,
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=1.0, By=0.0, Bz=3.0),
    )
    field_source = SumField(dipole, external)

    @test @inferred(field_source(0.0, 0.0, 0.0, 0.0)) ==
          EMField(SA[0.0, 0.0, 0.0], SA[1.0, 2.0, 3.0])
    @test SumField(()) isa ZeroField
    @test SumField((ZeroField(), dipole)) === dipole

    nested = SumField(ZeroField(), SumField(dipole, external), ZeroField())
    @test nested isa SumField
    @test nested.field_sources == (dipole, external)

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
    batch_field_source = MultipoleField(
      SA[1],
      SA[BatchParam([2.0, 3.0])],
      SA[BatchParam(0.0)],
    )
    lowered_batch_field_source = BeamTracking.batch_lower(batch_field_source)
    @test BeamTracking.static_batchcheck(lowered_batch_field_source)
    first_batch_field_source = @inferred BeamTracking.beval(lowered_batch_field_source, 1)
    second_batch_field_source = @inferred BeamTracking.beval(lowered_batch_field_source, 2)
    @test @inferred(first_batch_field_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
    @test @inferred(second_batch_field_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 3.0, 0.0]

    time_field_source = MultipoleField(SA[1], SA[4.0 * Time()], SA[TimeDependentParam(0.0)])
    evaluated_time_field_source =
      @inferred BeamTracking.teval(BeamTracking.time_lower(time_field_source), 0.5)
    @test @inferred(evaluated_time_field_source(0.0, 0.0, 0.0, 0.0)).B == SA[0.0, 2.0, 0.0]
  end

  @testset "No scalar allocations" begin
    field_source = SumField(
      MultipoleField(SA[2], SA[4.0], SA[5.0]),
      FunctionalField(
        test_uniform_field,
        (Ex=0.0, Ey=0.0, Ez=0.0, Bx=1.0, By=0.0, Bz=3.0),
      ),
    )
    @test_opt field_source(0.2, 0.3, 0.0, 0.0)
    @test @ballocated($field_source(0.2, 0.3, 0.0, 0.0)) == 0
  end

  @testset "Field source preparation" begin
    ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
    context = Beamlines.Context(strength=0.25)
    parameters = (
      strength=Beamlines.DefExpr{Float64}(c -> 2 * c.strength),
      constants=(label="map", count=2),
    )
    field_source = FunctionalField((x, y, s, t, p) -> begin
      v = zero(x)
      return EMField(v, v, v, v, v + p.strength, v)
    end, parameters)
    prepared = @inferred Beamlines.deval(field_source, context)
    @test prepared.evaluator === field_source.evaluator
    @test prepared.parameters.strength == 0.5
    @test prepared.parameters.constants === parameters.constants
    @test field_source.parameters.strength isa Beamlines.DefExpr

    opaque_parameters = FieldTestParameters(parameters.strength)
    opaque_field_source = FunctionalField(field_source.evaluator, opaque_parameters)
    @test Beamlines.deval(opaque_field_source, context).parameters === opaque_parameters
    @test BeamTracking.batch_lower(
      FunctionalField(field_source.evaluator, FieldTestParameters(BatchParam([0.5, 1.5]))),
    ).parameters isa FieldTestParameters{BatchParam}

    context.strength = 0.75
    @test Beamlines.deval(field_source, context).parameters.strength == 1.5

    dual_field_source = FunctionalField(field_source.evaluator, (
      strength=ForwardDiff.Dual(0.5, 1.0),
      constants=parameters.constants,
    ))
    scalar_field_source = @inferred Beamlines.scalarize(dual_field_source)
    @test scalar_field_source.evaluator === field_source.evaluator
    @test scalar_field_source.parameters.strength === 0.5

    dynamic = FunctionalField(field_source.evaluator, (
      strength=BatchParam([0.5, 1.5]),
      constants=parameters.constants,
    ))
    lowered = BeamTracking.batch_lower(dynamic)
    @test @inferred(BeamTracking.static_batchcheck(lowered))
    selected = @inferred BeamTracking.beval(lowered, 2)
    @test selected(0.0, 0.0, 0.0, 0.0).B[2] == 1.5
    @test_opt BeamTracking.beval(lowered, 2)
    @test @ballocated(BeamTracking.beval($lowered, 2)) == 0

    timed = FunctionalField(field_source.evaluator, (
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
    @test Beamlines.deval(opaque, context) === opaque
  end

  @testset "Numeric lowering" begin
    multipole = MultipoleField(SA[1, 2], SA[0.01, 0.03], SA[0.0, 0.02])
    functional = FunctionalField(test_uniform_field, (
      Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=0.1, Bz=0.0,
      nested=(values=(0.25, SA[0.5, 0.75]), order=2, label="field"),
    ))
    field_source = SumField(multipole, functional)
    for T in (Float32, Float16)
      lowered = @inferred BeamTracking.num_lower(T, field_source)
      @test lowered.field_sources[1].orders === multipole.orders
      @test lowered.field_sources[2].evaluator === functional.evaluator
      parameters = lowered.field_sources[2].parameters
      @test parameters.nested.values === (T(0.25), SVector{2,T}(0.5, 0.75))
      @test parameters.nested.order === 2
      @test parameters.nested.label === "field"
      @test @inferred(lowered(zero(T), zero(T), zero(T), zero(T))) isa EMField{T}
    end
    @test BeamTracking.num_lower(Float64, field_source) === field_source
    @test functional.parameters.By === 0.1
    opaque = FunctionalField(test_uniform_field, FieldTestParameters(0.1))
    @test BeamTracking.num_lower(Float32, opaque).parameters === opaque.parameters

    # Exercise the actual order: make_kernel_call lowers batch/time wrappers,
    # then pushing into a chain converts numbers to the coordinate precision.
    batch = MultipoleField(SA[1], SA[BatchParam([0.1, 0.2])], SA[BatchParam(0.0)])
    timed = FunctionalField(test_uniform_field, (
      Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0,
      By=TimeDependentParam(t -> Float32(t), false), Bz=0.0,
    ))
    call = BeamTracking.make_kernel_call(identity, (SumField(batch, timed),))
    chain = BeamTracking.KernelChain(Val{1}(), BeamTracking.RefState{Float32}(;
      t_enter=0f0, beta_gamma_enter=1f0))
    prepared = BeamTracking.push(chain, call).chain[1].args[1]
    for i in 1:2
      evaluated = BeamTracking.teval(BeamTracking.beval(prepared, i), 0.25f0)
      field = @inferred evaluated(0f0, 0f0, 0f0, 0f0)
      @test field isa EMField{Float32}
      @test field.B[2] ≈ Float32(i * 0.1) + 0.25f0
    end
    @test batch.normal[1].batch == [0.1, 0.2]
    bad_time = FunctionalField(test_uniform_field, (By=Time(),))
    bad_call = BeamTracking.make_kernel_call(identity, (bad_time,))
    @test_throws ErrorException BeamTracking.push(chain, bad_call)
  end

  @testset "Adaptation" begin
    multipole = MultipoleField(SA[0, 2], SA[0.01, 0.03], SA[0.0, 0.02])
    adapted_multipole = BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), multipole)
    @test adapted_multipole == multipole
    @test @ballocated(BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), $multipole)) == 0

    field_source = SumField(
      MultipoleField(SA[1], SA[0.01], SA[0.0]),
      FunctionalField(test_uniform_field, (field_map=[1.0, 2.0, 3.0],)),
    )
    adapted = BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), field_source)

    @test adapted isa SumField
    @test adapted.field_sources[1] isa MultipoleField
    @test adapted.field_sources[2] isa FunctionalField
    @test adapted.field_sources[2].parameters.field_map == SA[1.0, 2.0, 3.0]
  end
end

@testset "Field unit conventions" begin
  for T in (Float32, Float64), R in (T(-3), T(2))
    args = (T(0.02), T(-0.01), zero(T), zero(T))
    physical = FunctionalField((x,y,s,t,p) -> EMField(p...), T.((1,2,3,4,5,6)))
    normalized = FunctionalField(physical.evaluator, physical.parameters ./ R; normalized=true)
    expected = physical(args...)
    converted = @inferred BeamTracking.normalized_field_at(physical, args..., inv(R))
    direct = @inferred BeamTracking.normalized_field_at(normalized, args..., inv(R))
    @test converted.E ≈ expected.E / R
    @test converted.B ≈ expected.B / R
    @test converted.E ≈ direct.E
    @test converted.B ≈ direct.B
    @test physical(args...).E == expected.E
    parameter_free = FunctionalField((x,y,s,t) -> EMField(x,x,x,y,y,y); normalized=true)
    @test @inferred(BeamTracking.normalized_field_at(parameter_free, args..., inv(R))) == parameter_free(args...)
    multipole = MultipoleField(SA[1,2], T.(SA[0.1,0.2]), T.(SA[0,0]); normalized=true)
    mixed = SumField(multipole, physical)
    result = @inferred BeamTracking.normalized_field_at(mixed, args..., inv(R))
    @test result.E ≈ converted.E
    @test result.B ≈ multipole(args...).B + converted.B
    @test_throws ArgumentError mixed(args...)
    same_units = SumField(multipole, normalized)
    @test same_units(args...).B ≈ result.B
    for field_source in (multipole, normalized, parameter_free)
      @test BeamTracking.field_normalized(BeamTracking.num_lower(Float32, field_source)) == Val(true)
      @test BeamTracking.field_normalized(BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), field_source)) == Val(true)
      @test BeamTracking.field_normalized(BeamTracking.time_lower(field_source)) == Val(true)
      @test BeamTracking.field_normalized(BeamTracking.batch_lower(field_source)) == Val(true)
    end
  end
end
