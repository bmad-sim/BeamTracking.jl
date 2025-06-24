
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
    ele_sol = LineElement(L=1.0, Ks=2.0, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele_sol], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)


    # Errors:
    ele_kick = LineElement(L=1.0, K0L=1.0, tracking_method=Exact())
    ele_bend = LineElement(L=1.0, g=0.01, tracking_method=Exact())
    ele_patch_bend = LineElement(L=1.0, g=0.01, dy=3.0, dz_rot=0.3, tracking_method=Exact())
    ele_patch_sol = LineElement(L=1.0, Ks=1.0, dt=1.0, tracking_method=Exact())
    ele_bend_quad = LineElement(L=1.0, g=0.01, K1=1.0, tracking_method=Exact())
    @test_throws ErrorException track!(b0, Beamline([ele_kick],       Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend],       Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_bend], Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_sol],  Brho_ref=Brho_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend_quad],  Brho_ref=Brho_ref))
  end

  @testset "Integration" begin
    p0c = 15e9
    # E to Brho
    Brho_ref = BeamTracking.calc_Brho(ELECTRON, sqrt(p0c^2 + BeamTracking.massof(ELECTRON)^2))

    # Thin pure multipole:
    ele = LineElement(L=0.0, K3L=10.0, tilt3=0.5*pi, tracking_method=Standard())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_pure_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Thin multipole:
    ele = LineElement(L=0.0, K2L=1.0, tilt2=0.3*pi, K6L=100.0, tilt6=0.15*pi, tracking_method=Standard())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/thin_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    p0c = 10e6
    # E to Brho
    Brho_ref = BeamTracking.calc_Brho(ELECTRON, sqrt(p0c^2 + BeamTracking.massof(ELECTRON)^2))

    # Drift: 
    ele = LineElement(L=1.0, tracking_method=Standard())   
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/drift.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure solenoid:
    ele = LineElement(L=1.0, Ks=2.0, tracking_method=Standard())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Solenoid with quadrupole:
    ele = LineElement(L=2.0, Ks=0.1, K1=0.1, tracking_method=Standard())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sol_quad.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-5)

    # Pure dipole:
    #ele = LineElement(L=2.0, K0=0.1+1e-3, g=0.1, e1=1e-3, e2=-2e-3, tracking_method=Standard())
    #b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    #bl = Beamline([ele], Brho_ref=Brho_ref)
    #track!(b0, bl)
    #v_expected = read_map("bmad_maps/sbend.jl")
    #@test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Dipole with quadrupole:
    #ele = LineElement(L=2.0, K0=0.1+1e-3, K1=0.2, g=0.1, e1=1e-3, e2=-2e-3, tracking_method=Standard())
    #b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    #bl = Beamline([ele], Brho_ref=Brho_ref)
    #track!(b0, bl)
    #v_expected = read_map("bmad_maps/sbend_quad.jl")
    #@test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure quadrupole:
    ele = LineElement(L=2.0, K1=0.1, tracking_method=Standard(order=6,num_steps=100))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-5)

    # Skewed quadrupole:
    ele = LineElement(L=2.0, K1=0.1, tilt1 = 0.25*pi, tracking_method=Standard(order=6,num_steps=100))
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/skewed_quad.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-5)

    # Quadrupole with octupole:
    ele = LineElement(L=2.0, K1=0.1, K3=100, tracking_method=DKD())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/quad_oct.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 5e-10)

    # Pure sextupole:
    ele = LineElement(L=2.0, K2=0.1, tracking_method=Standard())
    b0 = Bunch(collect(transpose(@vars(D10))), Brho_ref=Brho_ref)
    bl = Beamline([ele], Brho_ref=Brho_ref)
    track!(b0, bl)
    v_expected = read_map("bmad_maps/sextupole.jl")
    @test coeffs_approx_equal(v_expected, b0.v, 1e-5)
  end
end