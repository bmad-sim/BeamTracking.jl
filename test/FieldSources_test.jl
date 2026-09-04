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

struct FieldTestStrength{T}
  value::T
end

struct FieldTestParameters{S,G}
  strength::S
  grid::G
end

struct FieldTestCustomSource{P}
  parameters::P
end

BeamTracking.Adapt.@adapt_structure FieldTestStrength
BeamTracking.Adapt.@adapt_structure FieldTestParameters
BeamTracking.Adapt.@adapt_structure FieldTestCustomSource

struct FieldTestOpaque{T}
  value::T
end
BeamTracking.field_parameter_leaf(::Type{<:FieldTestOpaque}) = true

function (source::FieldTestCustomSource)(x, y, z, s)
  v = zero(x)
  return EMField(v, v, v, v, v + source.parameters.strength.value, v)
end

# The dimension is not inferable from fields, and construction is keyword-only.
struct FieldTestSpecialSource{T,N}
  strength::T
  function FieldTestSpecialSource(; strength::T, dimension) where {T}
    return new{T,dimension}(strength)
  end
end

function (source::FieldTestSpecialSource)(x, y, z, s)
  v = zero(x)
  return EMField(v, v, v, v, v + source.strength, v)
end

BeamTracking.rebuild_field_source(::FieldTestSpecialSource{T,N}, children::Tuple) where {T,N} =
  FieldTestSpecialSource(strength=only(children), dimension=N)

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

  @testset "Recursive source preparation" begin
    ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
    context = Beamlines.Context(strength=0.25)
    grid = [1.0, 2.0, 3.0]
    parameters = FieldTestParameters(
      FieldTestStrength(Beamlines.DefExpr{Float64}(c -> 2 * c.strength)), grid,
    )
    source = FieldTestCustomSource(parameters)
    prepared = @inferred ext.unpack_field_source(source, context)
    @test prepared.parameters.strength isa FieldTestStrength{Float64}
    @test prepared.parameters.strength.value == 0.5
    @test prepared.parameters.grid === grid
    @test source.parameters.strength.value isa Beamlines.DefExpr

    opaque = FieldTestOpaque(parameters.strength.value)
    @test @inferred(ext.unpack_field_source(opaque, context)) === opaque

    context.strength = 0.75
    @test ext.unpack_field_source(source, context).parameters.strength.value == 1.5

    special = FieldTestSpecialSource(strength=parameters.strength.value, dimension=3)
    prepared_special = @inferred ext.unpack_field_source(special, context)
    @test prepared_special isa FieldTestSpecialSource{Float64,3}
    @test prepared_special.strength == 1.5

    dual_source = FieldTestCustomSource(FieldTestParameters(
      FieldTestStrength(ForwardDiff.Dual(0.5, 1.0)), grid,
    ))
    scalar_source = @inferred ext.scalarize_field_source(dual_source)
    @test scalar_source.parameters.strength.value === 0.5

    # FunctionalField preserves even stateful evaluator objects unchanged.
    functional = FunctionalField(special, parameters)
    prepared_functional = @inferred ext.unpack_field_source(functional, context)
    @test prepared_functional.evaluator === special
    @test prepared_functional.parameters.strength.value == 1.5

    closure = let deferred = parameters.strength.value
      () -> deferred
    end
    @test ext.unpack_field_source(closure, context) === closure
    @test isequal(ext.unpack_field_source((label="map", count=2, missing=missing), context),
                  (label="map", count=2, missing=missing))

    # Use the same generated traversal for nested batch and time parameters.
    dynamic = FieldTestCustomSource(FieldTestParameters(
      FieldTestStrength(BatchParam([0.5, 1.5])), grid,
    ))
    wrapped = BeamTracking._PreparedField(dynamic)
    lowered = BeamTracking.batch_lower(wrapped)
    @test @inferred(BeamTracking.static_batchcheck(lowered))
    selected = @inferred BeamTracking.beval(lowered, 2)
    @test selected(0.0, 0.0, 0.0, 0.0).B[2] == 1.5
    @test selected.source.parameters.grid === grid
    @test_opt BeamTracking.beval(lowered, 2)
    @test @ballocated(BeamTracking.beval($lowered, 2)) == 0

    timed = BeamTracking._PreparedField(FieldTestCustomSource(FieldTestParameters(
      FieldTestStrength(2 * Time()), grid,
    )))
    lowered_time = BeamTracking.time_lower(timed)
    @test @inferred(BeamTracking.static_timecheck(lowered_time))
    evaluated = @inferred BeamTracking.teval(lowered_time, 0.25)
    @test evaluated(0.0, 0.0, 0.0, 0.0).B[2] == 0.5
    @test @ballocated(BeamTracking.teval($lowered_time, 0.25)) == 0

    # Nested SVector parameters must lower to tuples before SIMD evaluation.
    vector_source = FunctionalField(test_uniform_field, (
      coefficients=SA[FieldTestStrength(BatchParam([0.5, 1.5]))],
      constants=SA[1.0, 2.0],
    ))
    vector_lowered = BeamTracking.batch_lower(vector_source)
    @test vector_lowered.parameters.coefficients isa Tuple
    @test vector_lowered.parameters.constants === SA[1.0, 2.0]
    vector_selected = @inferred BeamTracking.beval(vector_lowered, SIMD.VecRange{2}(1))
    @test Tuple(vector_selected.parameters.coefficients[1].value) == (0.5, 1.5)

    time_vector = FunctionalField(test_uniform_field, (coefficients=SA[Time(), 2 * Time()],))
    time_vector_lowered = BeamTracking.time_lower(time_vector)
    @test time_vector_lowered.parameters.coefficients isa Tuple
    @test @inferred(BeamTracking.teval(time_vector_lowered, 0.25)).parameters.coefficients ==
          (0.25, 0.5)

    adapted = BeamTracking.Adapt.adapt(FieldSourceTestAdaptor(), BeamTracking._PreparedField(prepared))
    @test adapted.source.parameters.grid == SA[1.0, 2.0, 3.0]
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
