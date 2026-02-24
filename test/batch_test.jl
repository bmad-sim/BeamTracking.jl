using Beamlines
using BeamTracking
using Test
using LinearAlgebra

# Pushes 4 equal particles with 4 batched parameters
# compares when particle is pushed through 4 different elements
function test_batch(
  batch_params, 
  nonbatch_params=(;); 
  v=nothing, 
  q=nothing, 
  extra_tests=(),
)

  if isnothing(v)
    v = repeat((rand(1,6).-0.5)*1e-4, 4, 1) 
  end
  if isnothing(q)
    q = repeat(mapslices(normalize, rand(1,4), dims=2), 4, 1)
  end

  E_ref = 5e9
  species = Species("electron")

  ele_1 = LineElement()
  ele_2 = LineElement()
  ele_3 = LineElement()
  ele_4 = LineElement()
  ele_batch = LineElement()
  for (k,val) in pairs(batch_params)
    val_1 = -val
    val_2 = -0.5*val
    val_3 = 0.5*val
    val_4 = val
    setproperty!(ele_1, k, val_1)
    setproperty!(ele_2, k, val_2)
    setproperty!(ele_3, k, val_3)
    setproperty!(ele_4, k, val_4)
    setproperty!(ele_batch, k, BatchParam([val_1, val_2, val_3, val_4]))
  end
  for (k,val) in pairs(nonbatch_params)
    setproperty!(ele_1, k, val)
    setproperty!(ele_2, k, val)
    setproperty!(ele_3, k, val)
    setproperty!(ele_4, k, val)
    setproperty!(ele_batch, k, val)
  end

  bl_1 = Beamline([ele_1], E_ref=E_ref, species_ref=species)
  bl_2 = Beamline([ele_2], E_ref=E_ref, species_ref=species)
  bl_3 = Beamline([ele_3], E_ref=E_ref, species_ref=species)
  bl_4 = Beamline([ele_4], E_ref=E_ref, species_ref=species)
  bl_batch = Beamline([ele_batch], E_ref=E_ref, species_ref=species)

  b0_1 = Bunch(v[1,:]', q[1,:]'; p_over_q_ref=bl_1.p_over_q_ref, species=bl_1.species_ref)
  b0_2 = Bunch(v[1,:]', q[1,:]'; p_over_q_ref=bl_2.p_over_q_ref, species=bl_2.species_ref)
  b0_3 = Bunch(v[1,:]', q[1,:]'; p_over_q_ref=bl_3.p_over_q_ref, species=bl_3.species_ref)
  b0_4 = Bunch(v[1,:]', q[1,:]'; p_over_q_ref=bl_4.p_over_q_ref, species=bl_4.species_ref)

  b0_batch = Bunch(v, q; p_over_q_ref=bl_batch.p_over_q_ref, species=bl_batch.species_ref)

  track!(b0_1, bl_1)
  track!(b0_2, bl_2)
  track!(b0_3, bl_3)
  track!(b0_4, bl_4)

  # Ensure branchlessness of parameters with explicit SIMD
  if (VERSION < v"1.11" && Sys.ARCH == :x86_64)
    use_explicit_SIMD=false
  else
    use_explicit_SIMD=true
  end
  track!(b0_batch, bl_batch; use_explicit_SIMD=use_explicit_SIMD) 
  
  @test b0_batch.coords.v[1,:]' ≈ b0_1.coords.v
  @test b0_batch.coords.v[2,:]' ≈ b0_2.coords.v
  @test b0_batch.coords.v[3,:]' ≈ b0_3.coords.v
  @test b0_batch.coords.v[4,:]' ≈ b0_4.coords.v

  @test b0_batch.coords.q[1,:]' ≈ b0_1.coords.q
  @test b0_batch.coords.q[2,:]' ≈ b0_2.coords.q
  @test b0_batch.coords.q[3,:]' ≈ b0_3.coords.q
  @test b0_batch.coords.q[4,:]' ≈ b0_4.coords.q
  
  for extra_test in extra_tests
    @test extra_test(b0_plus, b0_minus, b0_time)
  end
end

@testset "Batch" begin
  # Test each of the splits: DKD, MKM, SKS, BKB
  # MKM:
  test_batch((;Kn1=0.36), (;L=0.5))
  test_batch((;Kn1=0.36, Ks2=-1.2, Kn12=105.), (;L=3.4))

  # SKS:
  test_batch((;Ksol=1.2), (;L=3.4)) 
  test_batch((;Ksol=0.23, Kn1=0.36, Ks2=-1.2, Kn12=105.), (;L=3.4))

  # DKD:
  test_batch((;Kn0=1e-2, Kn1=-2, Kn3=14), (;L=5.1))

  # BKB not working at the moment because bend multipoles not 
  # implemented and Ks0 (also stored) becomes a TimeDependentParam
  # test_batch((;Kn0=1e-4, e1=2e-2, e2=3e-2), (;L=1.4))
  
  Kn1 = 0.36
  Ks2 = -1.2
  L = 3.4
  Kn12 = 105.
  Kn3 = 14
  Kn0 = 1e-2
  Ksol = -0.23
  test_batch((;Kn1=Kn1), (;L=L))
  test_batch((;Kn1=Kn1, Ks2=Ks2, Kn12=Kn12), (;L=L))

  # SKS:
  test_batch((;Ksol=Ksol), (;L=L))
  test_batch((;Ksol=Ksol, Kn1=Kn1, Ks2=Ks2, Kn12=Kn12), (;L=L))

  # DKD:
  test_batch((;Kn0=Kn0, Kn1=Kn1, Kn3=Kn3), (;L=L))

  # Now with different types of multipoles entered:
  test_batch((;Bn1=Kn1), (;L=0.5))
  test_batch((;Bn1L=Kn1, Ks2L=Ks2, Bn12=Kn12), (;L=L))

  # SKS:
  test_batch((;Bsol=Ksol), (;L=L))
  test_batch((;BsolL=Ksol), (;L=L))
  test_batch((;KsolL=Ksol, Bn1L=Kn1, Bs2=Ks2, Kn12L=Kn12), (;L=L))

  # DKD:
  test_batch((;Bn0L=Kn0, Kn1L=Kn1, Bn3=Kn3), (;L=L))
  test_batch((;Bn0L=Kn0, Kn1L=Kn1, Bn3=Kn3), (;L=L))
  test_batch((;Bn0L=Kn0, Kn1L=Kn1, Bn3=Kn3), (;L=0))
  test_batch((;Bn0L=Kn0, Kn1L=Kn1, Bn3=Kn3), (;L=0))

  #=
  # Aperture:
  # let's make a time-dependent aperture which oscillates but will allow both
  # particles through
  function check_state(p, m, t, state=BeamTracking.STATE_ALIVE)
    if !all(p.coords.state .== m.coords.state .== t.coords.state .== state)
      error("Test failed: p.coords.state = $(p.coords.state), m.coords.state = $(m.coords.state),
      t.coords.state = $(t.coords.state), state = $state")
    else
      return true
    end
  end

  test_batch(
    (;x2_limit=1), 
    (;aperture_shape = ApertureShape.Rectangular); 
    v=[0.5 0 0 0 0 0],
    ft=(v)->v*cos(2*Time()),
    extra_tests=(check_state,)
  )

  # now one where all particle should die
  test_batch(
    (;x2_limit=1), 
    (;aperture_shape = ApertureShape.Rectangular); 
    v=[0.5 0 0 0 0 0],
    ft=(v)->v*sin(2*Time()),
    extra_tests=((p,m,t)->check_state(p,m,t,BeamTracking.STATE_LOST_POS_X),)
  )

=#

  # Test that batch evaluation works for structs
  # This tests batch_lower, beval, and static_batchcheck on structs with BatchParam fields
  @testset "Batch evaluation for structs" begin
    # Create a simple struct with flexible types - one BatchParam field, one constant
    struct TestStruct{T1,T2}
      x::T1
      y::T2
    end

    # Create test data
    batch_values_x = [1.0, 2.0, 3.0, 4.0]
    batch_x = BatchParam(batch_values_x)
    test_struct = TestStruct(batch_x, 5.0)
    
    # Test 1: batch_lower converts BatchParam to _LoweredBatchParam
    # constant field remains unchanged
    lowered = BeamTracking.batch_lower(test_struct)
    @test lowered.x isa BeamTracking._LoweredBatchParam
    @test lowered.y == 5.0
    @test length(lowered.x.batch) == 4
    
    # Test 2: static_batchcheck detects _LoweredBatchParam in lowered struct
    @test BeamTracking.static_batchcheck(lowered) == true
    
    # Test 3: beval evaluates struct per particle index
    for i in 1:4
      evaluated = BeamTracking.beval(lowered, i)
      @test evaluated isa TestStruct
      @test evaluated.x == batch_values_x[i]
      @test evaluated.y == 5.0
    end
    
    # Test 4: No batch params - static_batchcheck returns false
    test_struct_no_batch = TestStruct(1.0, 2.0)
    @test BeamTracking.static_batchcheck(test_struct_no_batch) == false

    # Test 5: batch_lower returns same object when no fields changed
    @test BeamTracking.batch_lower(test_struct_no_batch) === test_struct_no_batch
    
    # Test 5: Nested structs
    struct OuterStruct{T}
      inner::T
      scale::Float64
    end
    
    nested_struct = OuterStruct(test_struct, 2.0)
    lowered_nested = BeamTracking.batch_lower(nested_struct)
    @test lowered_nested.inner.x isa BeamTracking._LoweredBatchParam
    @test BeamTracking.static_batchcheck(lowered_nested) == true
    
    for i in 1:4
      evaluated = BeamTracking.beval(lowered_nested, i)
      @test evaluated isa OuterStruct
      @test evaluated.inner isa TestStruct
      @test evaluated.inner.x == batch_values_x[i]
      @test evaluated.inner.y == 5.0
      @test evaluated.scale == 2.0
    end
  end
end