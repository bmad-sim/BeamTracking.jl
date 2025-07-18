@testset "Beamlines" begin
  include("lattices/esr.jl")


  @testset "Linear" begin
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=ring.Brho_ref)
    foreach(t->t.tracking_method=Linear(), ring.line)
    track!(b0, ring)
    M_ESR = [  0.8763088913632391E+00  0.2842332738570903E+00 -0.9233408598814828E-06 -0.1104742931103878E-06  0.0000000000000000E+00 -0.8843595261589663E-07
              -0.8165279324836310E+00  0.8763069736287854E+00 -0.1898521129265218E-05 -0.1113630193834745E-06  0.0000000000000000E+00 -0.1417461685411299E-06
               0.1352460571822227E-07 -0.2969583000777938E-07  0.6374265424413608E+00  0.4391919124687460E-01  0.0000000000000000E+00  0.5592919170820314E-14
              -0.4079601917518316E-05  0.1052556482124766E-05 -0.1351773220872402E+02  0.6374258159142916E+00  0.0000000000000000E+00 -0.4735132312703365E-13
               0.1964238541540386E-06 -0.3720806396318763E-07 -0.8403098936929155E-14 -0.1661103086541925E-15  0.1000000000000000E+01 -0.2358669003370585E+01
               0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  ]

    @test GTPSA.jacobian(b0.v) ≈ M_ESR

    bblring = BitsBeamline(ring, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=ring.Brho_ref)
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.v) ≈ M_ESR

    #combined
    ring2 = Beamline([SBend(Kn0 = 0.24, g = 0.25, Kn1 = 0.25, L = 1.0, e1 = 0.5, e2 = 0.2)], Brho_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=-0.0017045090263411496)
    foreach(t -> t.tracking_method = Linear(), ring2.line)
    track!(b0, ring2)
    M_combined = [
      0.9734056928978285E+00  0.9491282811321522E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1170583759019584E+00
      -0.1355634359386695E+00  0.8951384971554924E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.2434778344045677E+00
       0.0000000000000000E+00  0.0000000000000000E+00  0.9909816470445496E+00  0.1042190610987494E+01  0.0000000000000000E+00  0.0000000000000000E+00
       0.0000000000000000E+00  0.0000000000000000E+00  0.6449002140313032E-01  0.1076922966223960E+01  0.0000000000000000E+00  0.0000000000000000E+00
      -0.2528715457465107E+00 -0.1263082397778435E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  0.4893946839414888E+00
       0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01
    ]
    @test GTPSA.jacobian(b0.v) ≈  M_combined

    bblring = BitsBeamline(ring2, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=ring2.Brho_ref)
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.v) ≈ M_combined

    #g=0, dg!=0 thin corrector coil

    #SBend(Kn0=0.17, g=0, L = 1.1), line #L72 - L76
    ring3 = Beamline([SBend(Kn0=0.17,g=0,Kn1=-0.20,L=1.1)], Brho_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=-0.0017045090263411496)
    foreach(t -> t.tracking_method = Linear(), ring3.line)
    track!(b0, ring3)
    M_combined = [
      0.1123459935969970E+01  0.1144906606954581E+01  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1070487677502533E+00
      0.2289813213909162E+00  0.1123459935969970E+01  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.7726442422052786E-02
      0.0000000000000000E+00  0.0000000000000000E+00  0.8814205674902971E+00  0.1056167087172011E+01  0.0000000000000000E+00  0.0000000000000000E+00
      0.0000000000000000E+00  0.0000000000000000E+00 -0.2112334174344023E+00  0.8814205674902971E+00  0.0000000000000000E+00  0.0000000000000000E+00
      0.1583181978396722E-01  0.1114189467851014E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  0.5637819347138261E+00
      0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01
    ]
    @test GTPSA.jacobian(b0.v) ≈  M_combined

    bblring = BitsBeamline(ring3, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=ring2.Brho_ref)
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.v) ≈ M_combined


    #thick_bend_no_field line 103-104 --> singularity
    #L132 - L134
    ring4 = Beamline([SBend(Kn0=0.0, g=0.005, L = 1.1, Kn1 = 0.05, e1 = 0.1, e2 = 0.2)], Brho_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=-0.0017045090263411496)
    foreach(t -> t.tracking_method = Linear(), ring4.line)
    track!(b0, ring4)
    M_combined = [
      0.9699022031860474E+00  0.1088941837334312E+01  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1518962872590519E-04
      -0.5444709186671562E-01  0.9699022031860474E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.5499832721245711E-02
       0.0000000000000000E+00  0.0000000000000000E+00  0.1030402818311841E+01  0.1111125267330245E+01  0.0000000000000000E+00  0.0000000000000000E+00
       0.0000000000000000E+00  0.0000000000000000E+00  0.5555626336651224E-01  0.1030402818311841E+01  0.0000000000000000E+00  0.0000000000000000E+00
      -0.5335126904601591E-02 -0.5974265494137833E-02  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  0.5499999184222805E+00
       0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01
    ]
    @test GTPSA.jacobian(b0.v) ≈  M_combined

    bblring = BitsBeamline(ring4, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), Brho_ref=ring2.Brho_ref)
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.v) ≈ M_combined

  end

  @testset "Exact" begin
    p0c = 10e6
    # E to Brho
    Brho_ref = BeamTracking.calc_Brho(ELECTRON, sqrt(p0c^2 + BeamTracking.massof(ELECTRON)^2))

    # Patch:
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele_patch], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Drift: 
    ele_drift = LineElement(L=1.0, tracking_method=Exact())   
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele_drift], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/drift.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)


    # Thick solenoid:
    ele_sol = LineElement(L=1.0, Ksol=2.0, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele_sol], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)


    # Errors:
    ele_kick = LineElement(L=1.0, Kn0L=1.0, tracking_method=Exact())
    ele_bend = LineElement(L=1.0, g=0.01, tracking_method=Exact())
    ele_patch_bend = LineElement(L=1.0, g=0.01, dy=3.0, dz_rot=0.3, tracking_method=Exact())
    ele_patch_sol = LineElement(L=1.0, Ksol=1.0, dt=1.0, tracking_method=Exact())
    ele_bend_quad = LineElement(L=1.0, g=0.01, Kn1=1.0, tracking_method=Exact())
    @test_throws ErrorException track!(b0, Beamline([ele_kick],       Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend],       Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_bend], Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_sol],  Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend_quad],  Brho_ref=Brho_ref))
  end

  @testset "SplitIntegration" begin
    p0c = 10e6
    # E to Brho
    Brho_ref = BeamTracking.calc_Brho(ELECTRON, sqrt(p0c^2 + BeamTracking.massof(ELECTRON)^2))

    # Thin straight pure dipole:
    ele = LineElement(L=0.0, Kn0L=1.0, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_pure_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Thin straight dipole:
    ele = LineElement(L=0.0, Kn0L=1.0, Kn5L=1000.0, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Thin pure quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1 = pi, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_pure_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Thin quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1 = pi, Kn5L=100.0, tilt5=-0.1*pi, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Thin pure multipole:
    ele = LineElement(L=0.0, Kn3L=10.0, tilt3=0.5*pi, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_pure_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Thin multipole:
    ele = LineElement(L=0.0, Kn2L=1.0, tilt2=0.3*pi, Kn6L=100.0, tilt6=0.15*pi, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-14)

    # Drift: 
    ele = LineElement(L=1.0, tracking_method=SplitIntegration())   
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/drift.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure solenoid:
    ele = LineElement(L=1.0, Ksol=2.0, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Solenoid with quadrupole:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sol_quad.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # SK multiple steps:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=4, num_steps=2))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Step size:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=4, ds_step=1.0))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Straight pure dipole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, tilt0=pi/3, tracking_method=DriftKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/straight_pure_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Straight dipole with quadrupole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.03, tracking_method=DriftKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/straight_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quadrupole_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Skew quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=pi/4, tracking_method=MatrixKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Skew quadrupole another way (MK):
    ele = LineElement(L=2.0, Ks1=-0.1, tracking_method=MatrixKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure quadrupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=0.1*pi, tracking_method=DriftKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quadrupole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Quadrupole with octupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quad_oct_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # MK multiple steps:
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(num_steps=2))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/mk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Quadrupole with octupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=DriftKick())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quad_oct_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure sextupole:
    ele = LineElement(L=2.0, Kn2=10.0, tilt2=0.2*pi, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sextupole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-9)

    # DK multiple steps:
    ele = LineElement(L=2.0, Kn2=10.0, tracking_method=SplitIntegration(order=4, num_steps=2))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/dk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Sextupole with decapole:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sex_dec.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Order 4:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=4))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/order_four.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Order 6:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=6))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/order_six.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-9)

    # Order 8:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=8))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/order_eight.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-9)

    # Patch:
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=SplitIntegration())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele_patch], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Errors:
    @test_throws ErrorException MatrixKick(ds_step = 0.1, num_steps = 2)
    @test_throws ErrorException BendKick(order = 2, num_steps = -2)
    @test_throws ErrorException DriftKick(ds_step = -0.1)
    @test_throws ErrorException SolenoidKick(num_steps = -2)
    @test_throws ErrorException SplitIntegration(order = 5)
  end

end