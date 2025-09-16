using Beamlines

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

  ele_plus = LineElement()
  ele_minus = LineElement()
  ele_time = LineElement()
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


  bl_plus = Beamline([ele_plus], E_ref=E_ref, species_ref=species)
  bl_minus = Beamline([ele_minus], E_ref=E_ref, species_ref=species)
  bl_time = Beamline([ele_time], E_ref=E_ref, species_ref=species)

  b0_plus = Bunch(deepcopy(v), deepcopy(q); R_ref=bl_plus.R_ref, species=bl_plus.species_ref)
  b0_minus = Bunch(deepcopy(v), deepcopy(q); R_ref=bl_minus.R_ref, species=bl_minus.species_ref)
  b0_time = Bunch(vcat(v,v), vcat(q,q); R_ref=bl_time.R_ref, species=bl_time.species_ref)
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
  
  # Linear tracking does not currently support TimeDependentParams

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


  
end