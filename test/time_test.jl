using Beamlines

# Unless the element already has one, the split it used with the old default tracking method
# (SciBmadStandard): solenoid -> SKS, quadrupole first -> MKM, dipole or higher first -> DKD,
# zero length -> Exact
function default_tracking_method(ele)
  isnothing(ele.tracking_method) || return ele.tracking_method
  ele.L == 0 && return Exact()
  isnothing(ele.BMultipoleParams) && return DriftKick()
  o = first(ele.BMultipoleParams.order)
  return o == 0 ? SolenoidKick() : o == 2 ? MatrixKick() : DriftKick()
end

# Given a list of element parameters, tracks
# two elements with same initial particle coordinates,
# but element parameters' signs flipped. Then compares
# that with tracking two particles through a single
# time dependent element where the parameters go with ft
# and the two particles' z coordinates are separated by C_LIGHT*pi
# so one particle sees positive strengths other negative
function test_time(
  dynamic_params, 
  static_params=(;); 
  v=nothing, 
  q=nothing, 
  ft=(v)->v*cos(Time()),
  extra_tests=(),
)

  if isnothing(v)
    v = (rand(1,6).-0.5)*1e-4; v[5] = 0; v[6] = 0;
  end
  if isnothing(q)
    q = rand(1,4)
    q = q/norm(q) # Normalize quaternion to 1
  end
  E_ref = 100e9
  species = Species("electron")

  ele_plus = LatElement()
  ele_minus = LatElement()
  ele_time = LatElement()
  for (k,val) in pairs(dynamic_params)
    setproperty!(ele_plus, k, Beamlines.deval(ft(val))(v[5]))
    setproperty!(ele_minus, k, Beamlines.deval(ft(val))(v[5]/BeamTracking.C_LIGHT-pi))
    setproperty!(ele_time, k, Beamlines.deval(ft(val)))
  end
  for (k,val) in pairs(static_params)
    setproperty!(ele_plus, k, val)
    setproperty!(ele_minus, k, val)
    setproperty!(ele_time, k, val)
  end


  ele_plus.tracking_method = default_tracking_method(ele_plus)
  ele_minus.tracking_method = default_tracking_method(ele_minus)
  ele_time.tracking_method = default_tracking_method(ele_time)

  bl_plus = Beamline([ele_plus], E_ref=E_ref, species_ref=species)
  bl_minus = Beamline([ele_minus], E_ref=E_ref, species_ref=species)
  bl_time = Beamline([ele_time], E_ref=E_ref, species_ref=species)

  b0_plus = Bunch(deepcopy(v), deepcopy(q); p_over_q_ref=bl_plus.p_over_q_ref, species=bl_plus.species_ref)
  b0_minus = Bunch(deepcopy(v), deepcopy(q); p_over_q_ref=bl_minus.p_over_q_ref, species=bl_minus.species_ref)
  b0_time = Bunch(vcat(v,v), vcat(q,q); p_over_q_ref=bl_time.p_over_q_ref, species=bl_time.species_ref)
  b0_time.coords.v[2,5] += -BeamTracking.C_LIGHT*pi

  track!(b0_plus, bl_plus)
  track!(b0_minus, bl_minus)
  # Ensure branchlessness of parameters with explicit SIMD
  track!(b0_time, bl_time; use_explicit_SIMD=true, use_KA=false) 
  
  @test b0_time.coords.v[1,1:4] ≈ b0_plus.coords.v[1:4]
  @test b0_time.coords.q[1,:]' ≈ b0_plus.coords.q
  @test b0_time.coords.v[2,1:4] ≈ b0_minus.coords.v[1:4]
  @test b0_time.coords.q[2,:]' ≈ b0_minus.coords.q

  for extra_test in extra_tests
    @test extra_test(b0_plus, b0_minus, b0_time)
  end
end

@testset "Time" begin
  # Right now length is not allowed to be a function of time 
  # It is kind of weird to do so but maybe possible so should 
  # be discussed.

  # Test each of the splits: DKD, MKM, SKS, BKB
  # MKM:
  test_time((;Kn1=0.36), (;L=0.5))
  test_time((;Kn1=0.36, Ks2=-1.2, Kn12=105.), (;L=3.4))

  # SKS:
  test_time((;Ksol=1.2), (;L=3.4)) 
  test_time((;Ksol=0.23, Kn1=0.36, Ks2=-1.2, Kn12=105.), (;L=3.4))

  # DKD:
  test_time((;Kn0=1e-2, Kn1=-2, Kn3=14), (;L=5.1))

  # BKB not working at the moment because bend multipoles not 
  # implemented and Ks0 (also stored) becomes a TimeDependentParam
  # test_time((;Kn0=1e-4, e1=2e-2, e2=3e-2), (;L=1.4))

  # Same but with DefExpr:
  Kn1 = 0.36
  Ks2 = -1.2
  L = 3.4
  Kn12 = 105.
  Kn3 = 14
  Kn0 = 1e-2
  Ksol = -0.23
  # MKM:
  test_time((;Kn1=DefExpr(()->Kn1)), (;L=DefExpr(()->L)))
  test_time((;Kn1=DefExpr(()->Kn1), Ks2=DefExpr(()->Ks2), Kn12=DefExpr(()->Kn12)), (;L=DefExpr(()->L)))

  # SKS:
  test_time((;Ksol=DefExpr(()->Ksol)), (;L=DefExpr(()->L)))
  test_time((;Ksol=DefExpr(()->Ksol), Kn1=DefExpr(()->Kn1), Ks2=DefExpr(()->Ks2), Kn12=DefExpr(()->Kn12)), (;L=DefExpr(()->L)))

  # DKD:
  test_time((;Kn0=DefExpr(()->Kn0), Kn1=DefExpr(()->Kn1), Kn3=DefExpr(()->Kn3)), (;L=DefExpr(()->L)))

  # Now with different types of multipoles entered:
  test_time((;Bn1=DefExpr(()->Kn1)), (;L=DefExpr(()->0.5)))
  test_time((;Bn1L=DefExpr(()->Kn1), Ks2L=DefExpr(()->Ks2), Bn12=DefExpr(()->Kn12)), (;L=DefExpr(()->3.4)))

  # SKS:
  test_time((;Bsol=DefExpr(()->Ksol)), (;L=DefExpr(()->3.4)))
  test_time((;BsolL=DefExpr(()->Ksol)), (;L=DefExpr(()->3.4)))
  test_time((;KsolL=DefExpr(()->Ksol), Bn1L=DefExpr(()->Kn1), Bs2=DefExpr(()->Ks2), Kn12L=DefExpr(()->Kn12)), (;L=DefExpr(()->3.4)))

  # DKD:
  test_time((;Bn0L=DefExpr(()->Kn0), Kn1L=DefExpr(()->Kn1), Bn3=DefExpr(()->Kn3)), (;L=DefExpr(()->L)))
  test_time((;Bn0L=Kn0, Kn1L=DefExpr(()->Kn1), Bn3=DefExpr(()->Kn3)), (;L=DefExpr(()->L)))
  test_time((;Bn0L=DefExpr(()->Kn0), Kn1L=DefExpr(()->Kn1), Bn3=DefExpr(()->Kn3)), (;L=0))
  test_time((;Bn0L=Kn0, Kn1L=DefExpr(()->Kn1), Bn3=DefExpr(()->Kn3)), (;L=0))

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

  test_time(
    (;x2_limit=1), 
    (;aperture_shape = ApertureShape.Rectangular); 
    v=[0.5 0 0 0 0 0],
    ft=(v)->v*cos(2*Time()),
    extra_tests=(check_state,)
  )

  # now one where all particle should die
  test_time(
    (;x2_limit=1), 
    (;aperture_shape = ApertureShape.Rectangular); 
    v=[0.5 0 0 0 0 0],
    ft=(v)->v*sin(2*Time()),
    extra_tests=((p,m,t)->check_state(p,m,t,BeamTracking.STATE_LOST_POS_X),)
  )

#=
  # Test that time evaluation works for structs
  # This tests time_lower, teval, and static_timecheck on structs with TimeDependentParam fields
  @testset "Time evaluation for structs" begin
    # Create a simple struct with flexible types - one TimeDependentParam field, one constant
    struct TimeTestStruct{T1,T2}
      x::T1
      y::T2
    end

    # Create test data
    tdp_x = 3.0 * Time()
    test_struct = TimeTestStruct(tdp_x, 5.0)

    # Test 1: time_lower converts TimeDependentParam to TimeFunction
    # constant field remains unchanged
    lowered = BeamTracking.time_lower(test_struct)
    @test lowered.x isa BeamTracking.TimeFunction
    @test lowered.y == 5.0

    # Test 2: static_timecheck detects TimeFunction in lowered struct
    @test BeamTracking.static_timecheck(lowered) == true

    # Test 3: teval evaluates struct at given time
    for t in [0.0, 1.0, 2.5, -1.0]
      evaluated = BeamTracking.teval(lowered, t)
      @test evaluated isa TimeTestStruct
      @test evaluated.x ≈ 3.0 * t
      @test evaluated.y == 5.0
    end

    # Test 4: No time params - static_timecheck returns false
    test_struct_no_time = TimeTestStruct(1.0, 2.0)
    @test BeamTracking.static_timecheck(test_struct_no_time) == false

    # Test 5: time_lower returns same object when no fields changed
    @test BeamTracking.time_lower(test_struct_no_time) === test_struct_no_time

    # Test 6: Nested structs
    struct TimeOuterStruct{T}
      inner::T
      scale::Float64
    end

    nested_struct = TimeOuterStruct(test_struct, 2.0)
    lowered_nested = BeamTracking.time_lower(nested_struct)
    @test lowered_nested.inner.x isa BeamTracking.TimeFunction
    @test BeamTracking.static_timecheck(lowered_nested) == true

    for t in [0.0, 1.0, 2.5]
      evaluated = BeamTracking.teval(lowered_nested, t)
      @test evaluated isa TimeOuterStruct
      @test evaluated.inner isa TimeTestStruct
      @test evaluated.inner.x ≈ 3.0 * t
      @test evaluated.inner.y == 5.0
      @test evaluated.scale == 2.0
    end
  end
=#
  
end
@testset "Float32 time dependence" begin
  # Float64 numbers captured by time functions are lowered with them, so Float32 tracking evaluates 
  # them in Float32 (also through literal_pow, whose closure has the static parameter N)
  tf = BeamTracking.time_lower(sqrt((0.3 + 1e5*Time())^2 + 0.2*Time()))
  tf32 = BeamTracking.num_lower(Float32, tf)
  @test tf32(1f-6) isa Float32
  @test tf32(1f-6) ≈ tf(1e-6) rtol=1e-6
  @test BeamTracking.num_lower(Float64, tf) === tf

  # Float32 tracking of a lattice with Float64 time dependent strengths (MatrixKick with a time 
  # dependent tilt, BendKick, DriftKick) and a per-particle ramp agrees with Float64 tracking
  e = Species("electron")
  R = Beamlines.E_to_R(e, 10e9)
  lat() = [Quadrupole(L=0.5, Kn1=-0.2 + 1e5*Time(), Ks1=0.05), Drift(L=1.0), 
           SBend(L=1.0, g=0.01, Kn1=0.01*(1 + 1e5*Time())), Sextupole(L=0.2, Kn2=1.0 + 1e5*Time())]
  v0 = 1e-3 .* [1 2 -1 0.5 3 1; -2 0 1 -1 0 -2; 0.5 -1 2 1 -3 0.5; 1 1 1 1 1 1]
  for (kw, p_over_q_ref) in (((;), R), ((; ramp_update_each_particle=true), R*(1 + 1e4*Time())))
    bl = Beamline(lat(); p_over_q_ref, species_ref=e)
    b64 = Bunch(copy(v0); species=e, p_over_q_ref=R)
    track!(b64, bl; kw...)
    for use_explicit_SIMD in (false, true)
      b32 = Bunch(Float32.(v0); species=e, p_over_q_ref=Float32(R))
      track!(b32, bl; use_explicit_SIMD, kw...)
      @test eltype(b32.coords.v) == Float32
      @test maximum(abs.(b32.coords.v .- b64.coords.v)) < 1e-6
      @test maximum(abs.(b64.coords.v .- v0)) > 1e-3 # the lattice does something
    end
  end
end
