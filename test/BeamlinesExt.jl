@testset "Beamlines" begin
  include("lattices/esr.jl")

  @testset "Linear" begin
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring.R_ref)
    foreach(t->t.tracking_method=Linear(), ring.line)
    track!(b0, ring)
    M_ESR = [  0.8763088913632391E+00  0.2842332738570903E+00 -0.9233408598814828E-06 -0.1104742931103878E-06  0.0000000000000000E+00 -0.8843595261589663E-07
              -0.8165279324836310E+00  0.8763069736287854E+00 -0.1898521129265218E-05 -0.1113630193834745E-06  0.0000000000000000E+00 -0.1417461685411299E-06
               0.1352460571822227E-07 -0.2969583000777938E-07  0.6374265424413608E+00  0.4391919124687460E-01  0.0000000000000000E+00  0.5592919170820314E-14
              -0.4079601917518316E-05  0.1052556482124766E-05 -0.1351773220872402E+02  0.6374258159142916E+00  0.0000000000000000E+00 -0.4735132312703365E-13
               0.1964238541540386E-06 -0.3720806396318763E-07 -0.8403098936929155E-14 -0.1661103086541925E-15  0.1000000000000000E+01 -0.2358669003370585E+01
               0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  ]

    @test GTPSA.jacobian(b0.coords.v) ≈ M_ESR

    bblring = BitsBeamline(ring, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring.R_ref, species=Species("electron"))
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.coords.v) ≈ M_ESR

    #combined
    ring2 = Beamline([SBend(Kn0 = 0.24, g_ref = 0.25, Kn1 = 0.25, L = 1.0, e1 = 0.5, e2 = 0.2)], R_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=-0.0017045090263411496, species=Species("electron"))
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
    @test GTPSA.jacobian(b0.coords.v) ≈  M_combined

    bblring = BitsBeamline(ring2, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring2.R_ref, species=Species("electron"))
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.coords.v) ≈ M_combined

    #g_ref=0, dg!=0 thin corrector coil

    #SBend(Kn0=0.17, g_ref=0, L = 1.1), line #L72 - L76
    ring3 = Beamline([SBend(Kn0=0.17,g_ref=0,Kn1=-0.20,L=1.1)], R_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=-0.0017045090263411496, species=Species("electron"))
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
    @test GTPSA.jacobian(b0.coords.v) ≈  M_combined

    bblring = BitsBeamline(ring3, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring2.R_ref, species=Species("electron"))
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.coords.v) ≈ M_combined


    #thick_bend_no_field line 103-104 --> singularity
    #L132 - L134
    ring4 = Beamline([SBend(Kn0=0.0, g_ref=0.005, L = 1.1, Kn1 = 0.05, e1 = 0.1, e2 = 0.2)], R_ref = -0.0017045090263411496)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=-0.0017045090263411496, species=Species("electron"))
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
    @test GTPSA.jacobian(b0.coords.v) ≈  M_combined

    bblring = BitsBeamline(ring4, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring2.R_ref, species=Species("electron"))
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.coords.v) ≈ M_combined

  end

  @testset "Exact" begin
    p0c = 10e6
    # E to R_ref
    R_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + BeamTracking.massof(Species("electron"))^2))
    
    # Patch:
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_patch], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Drift: 
    ele_drift = LineElement(L=1.0, tracking_method=Exact())   
    b0 = Bunch(collect(transpose(@vars(D10))), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_drift], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/drift.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Thick solenoid:
    ele_sol = LineElement(L=1.0, Ksol=2.0, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_sol], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Thick pure bend:
    ele_bend = LineElement(L=0.5, g_ref=1, Kn0 = 1.001, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    M_bend = [  0.8773527130168902E+00  0.4793669072377229E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1225248770305213E+00
               -0.4799049641428075E+00  0.8775825618903755E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.4794255386042034E+00
                0.0000000000000000E+00  0.0000000000000000E+00 0.1000000000000000E+01 0.4999794461108593E+00 0.0000000000000000E+00  0.0000000000000000E+00
                0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00
               -0.4794255937019158E+00 -0.1222950422117283E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 -0.1925165307287538E-01
                0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01   ]

    @test GTPSA.jacobian(b0.coords.v) ≈ M_bend

    # Thick pure dipole:
    ele_thick_dipole = LineElement(L=0.5, Kn0 = 0.001, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_thick_dipole], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    M_thick_dipole = [  0.1000000000000000E+01 0.5000000625000115E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1250000234375049E-03 
                        0.0000000000000000E+00 0.9999999999999988E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.1355252715606880E-18  
                        0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.5000000208333357E+00 0.0000000000000000E+00  0.0000000000000000E+00  
                        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00   
                        0.0000000000000000E+00 0.1250000234375048E-03 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01  0.1302241002743759E-02   
                        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01  ]

    @test GTPSA.jacobian(b0.coords.v) ≈ M_thick_dipole

    exact_bend_7 = 
      [ 0.9215672234619749E+00 0.4440044701830712E+00 0.0000000000000000E+00 -0.2757164812159411E-04 0.0000000000000000E+00 -0.9181358824490837E-01 
       -0.4314829847437838E+00 0.8772226326896770E+00 0.0000000000000000E+00 -0.1439716802793265E-03 0.0000000000000000E+00 -0.4794256953301575E+00  
        0.1445762786779847E-03 0.3973759587097068E-04 0.1000000000000000E+01  0.4547359879283919E+00 0.0000000000000000E+00  0.1302787444884017E-03  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00  
        0.4814390079976891E+00 0.1323261942503324E+00 0.0000000000000000E+00  0.1302787444884016E-03 0.1000000000000000E+01 -0.1960162611913873E-01 
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 ]
    ele_bend = LineElement(L=0.5, g_ref=-1, Kn0=-0.9, tracking_method=Exact())
    ps = [0.1, -7.5e-4, 0.1, -3e-4 , 0.1, -1e-3]
    b0 = Bunch(collect(transpose(@vars(D1) + ps)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_7

    # Large amplitude particles:
    exact_bend_8 = 
      [ 0.1080833710437112E+01  0.9094474670802266E+00 0.0000000000000000E+00  0.1550860894917999E+00 0.0000000000000000E+00 -0.3027871271030384E+00   
       -0.4799049641428072E+00  0.5214045791495372E+00 0.0000000000000000E+00 -0.3561779827408352E+00 0.0000000000000000E+00  0.6953951091606781E+00  
        0.3105425864451723E+00  0.4047877614568170E+00 0.1000000000000000E+01  0.1081052846155420E+01 0.0000000000000000E+00 -0.6712480744321767E+00 
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00  
       -0.6062974306786699E+00 -0.7902999152252140E+00 0.0000000000000000E+00 -0.6712480744321764E+00 0.1000000000000000E+01  0.5734407018255534E+00 
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 ]
    ele_bend = LineElement(L=0.5, g_ref=1, Kn0=1.001, tracking_method=Exact())
    ps2 = [0.9, 1.05, 0.9, 1.05, 0.9, 1.05]
    b0 = Bunch(collect(transpose(@vars(D1) + ps2)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_8

    exact_bend_9 = 
      [ -0.4593732760109923E+00 -0.4921986901333856E+00 0.0000000000000000E+00 -0.1930903588074987E+00 0.0000000000000000E+00 0.7723614352299949E+00  
        0.1513604990615857E+01 -0.5551163281719024E+00 0.0000000000000000E+00 0.1970545853834186E+00 0.0000000000000000E+00 -0.7882183415336744E+00 
        -0.2017409202902679E+00 0.2041776197971098E+00 0.1000000000000000E+01 0.2176998952198646E+01 0.0000000000000000E+00 0.5347268547649253E-01 
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00  
        0.8069636811610714E+00 -0.8167104791884392E+00 0.0000000000000000E+00 0.5347268547649256E-01 0.1000000000000000E+01 -0.2402983153080073E+01 
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ]
    ele_bend = LineElement(L=2, g_ref=2, Kn0=2, tracking_method=Exact())
    ps3 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps3)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_9

    exact_bend_10 = 
      [   0.1226742987625782E+01 0.1958192699228435E+01 0.0000000000000000E+00 0.6527308997428116E-01 0.0000000000000000E+00 -0.2610923598971246E+00 
          0.1875824874079783E-15 0.8151666731230962E+00 0.0000000000000000E+00 -0.1248317775345523E+00 0.0000000000000000E+00 0.4993271101382093E+00 
          0.1531365077233739E+00 0.2976531229986680E+00 0.1000000000000000E+01 0.1684214255208845E+01 0.0000000000000000E+00 -0.4582602041770446E+00 
          0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00
          -0.6125460308934955E+00 -0.1190612491994672E+01 0.0000000000000000E+00 -0.4582602041770448E+00 0.1000000000000000E+01 0.2646663249372617E+00 
          0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ]
    ele_bend = LineElement(L=2, g_ref=0.25, tracking_method=Exact())
    ps4 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps4)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_10

    exact_bend_11 = 
      [  0.9941038491254380E+00 0.1286686580424430E+01 0.0000000000000000E+00 0.4288955268081424E-01 0.0000000000000000E+00 -0.1715582107232570E+00 
         -0.2577684759557487E-16 0.1005931121662742E+01 0.0000000000000000E+00 0.5172908764299618E-01 0.0000000000000000E+00 -0.2069163505719847E+00  
         -0.5142408513764964E-01 -0.2341518705201763E-01 0.1000000000000000E+01 0.1356815341398884E+01 0.0000000000000000E+00 -0.3362769369682188E+00 
         0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00 
         0.2056963405505985E+00 0.9366074820807051E-01 0.0000000000000000E+00 -0.3362769369682188E+00 0.1000000000000000E+01 0.7363635310971241E-01 
         0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ]
    ele_bend = LineElement(L=2, g_ref=-0.1, tracking_method=Exact())
    ps5 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps5)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_11

    exact_bend_12 = 
      [ 0.1000000000000000E+01 0.1576112593915678E+01 0.0000000000000000E+00 0.1371117420537198E+00 0.0000000000000000E+00 -0.5484469682148793E+00  
        0.0000000000000000E+00 0.1000000000000005E+01 0.0000000000000000E+00 0.1016644466579460E-14 0.0000000000000000E+00 -0.4066577866317841E-14 
        0.0000000000000000E+00 0.1371117420537194E+00 0.1000000000000000E+01 0.1482335277082852E+01 0.0000000000000000E+00 -0.4202966917108471E+00 
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00 
        0.0000000000000000E+00 -0.5484469682148775E+00 0.0000000000000000E+00 -0.4202966917108473E+00 0.1000000000000000E+01 0.3052003750819157E+00 
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ]
    ele_bend = LineElement(L=2, Kn0=-0.3, tracking_method=Exact())
    ps6 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps6)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_12

    # With tilts:
    exact_bend_13 = 
      [  0.1121618635021878E+01 0.2324223731624664E+01 0.2106496550050672E+00 0.5324754046543224E+00 0.0000000000000000E+00 -0.4236468176166648E+00
         0.3767144394065392E-15 0.9195788093750110E+00 0.6524885489969529E-15 -0.1491420860171207E+00 0.0000000000000000E+00 0.2940407548192738E+00 
         0.2255432527764834E+00 0.5599237959506066E+00 0.1390652373113219E+01 0.3022460217047000E+01 0.0000000000000000E+00 -0.8256055287555776E+00 
         0.6524885489969529E-15 -0.1392935881676628E+00 0.1130143318219617E-14 0.7416783294715393E+00 0.0000000000000000E+00 0.5092935268428853E+00  
         -0.4446693087233498E+00 -0.1243157144846748E+01 -0.7701898352753723E+00 -0.2245039052845477E+01 0.1000000000000000E+01 0.7179291187108803E+00   
         0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ] 
    ele_bend = LineElement(L=2, g_ref=0.3, tilt_ref=pi/3, tracking_method=Exact())
    ps7 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps7)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_13

    exact_bend_14 = 
      [  0.1000000000000000E+01 0.1246718736520172E+01 -0.2561505167594759E-01 0.2187459575081172E-01 0.0000000000000000E+00 -0.1633116441631849E+00  
         -0.9683578852180224E-34 0.1000000000000000E+01 0.1581448440304960E-17 -0.1946922998101708E-17 0.0000000000000000E+00 0.1266997232096158E-16 
         0.3866162887417777E-17 0.5357085908743905E-01 0.9368607685071386E+00 0.1291327026000506E+01 0.0000000000000000E+00 -0.4285668726995124E+00   
         0.1581448440304961E-17 0.2586454382149800E-01 -0.2582701300335780E-01 0.1031795665484240E+01 0.0000000000000000E+00 -0.2069163505719838E+00   
         -0.1254775641797731E-16 -0.1633116441631849E+00 0.2049204134075807E+00 -0.1749967660064938E+00 0.1000000000000000E+01 0.8146308469937176E-01   
         0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ] 
    ele_bend = LineElement(L=2, g_ref=0.1, Kn0=0.13, tilt_ref=-pi/2, tracking_method=Exact())
    ps8 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps8)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_14

    # With skew strength:
    exact_bend_15 = 
      [ 0.1000000000000000E+01 0.1362928135021343E+01 0.0000000000000000E+00 0.6415149876723009E-01 0.0000000000000000E+00 -0.1922017848558854E+00 
        0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.1079213939202989E-30  
        0.0000000000000000E+00 0.6415149876723000E-01 0.1000000000000000E+01 0.1513589055740098E+01 0.0000000000000000E+00 -0.5132119901378400E+00  
        0.0000000000000000E+00 -0.2369523054950603E-15 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.1762490115442885E-14  
        0.0000000000000000E+00 -0.1922017848558854E+00 0.0000000000000000E+00 -0.5132119901378409E+00 0.1000000000000000E+01 0.1999860793263922E+00  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 ] 
    ele_bend = LineElement(L=2, Ks0=0.13, tracking_method=Exact())
    ps9 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    b0 = Bunch(collect(transpose(@vars(D1) + ps9)), R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_15

    # Particle lost (does not intersect exit face):
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_kick = LineElement(L=1.0, Kn0L=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_kick], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v

    # Particle lost (abs(px) > pt):
    b0 = Bunch([0.1 0.2 0.3 0.4 0.5 -0.6], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v

    # Particle lost in bend (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v

    # Particle lost in drift (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_drift = LineElement(L=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_drift], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v

    # Particle lost in solenoid (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_sol = LineElement(L=1.0, Ksol=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_sol], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v

    # Errors:
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_patch_bend = LineElement(L=1.0, g_ref=0.01, dy=3.0, dz_rot=0.3, tracking_method=Exact())
    ele_patch_sol = LineElement(L=1.0, Ksol=1.0, dt=1.0, tracking_method=Exact())
    ele_bend_quad = LineElement(L=1.0, g_ref=0.01, Kn1=1.0, tracking_method=Exact())
    @test_throws ErrorException track!(b0, Beamline([ele_patch_bend], R_ref=R_ref, species_ref=Species("electron")))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_sol],  R_ref=R_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend_quad],  R_ref=R_ref))
  end

  @testset "SplitIntegration" begin
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring.R_ref)
    foreach(t->t.tracking_method=SplitIntegration(), ring.line)
    track!(b0, ring)
    M_ESR = [  0.8763088913632391E+00  0.2842332738570903E+00 -0.9233408598814828E-06 -0.1104742931103878E-06  0.0000000000000000E+00 -0.8843595261589663E-07
              -0.8165279324836310E+00  0.8763069736287854E+00 -0.1898521129265218E-05 -0.1113630193834745E-06  0.0000000000000000E+00 -0.1417461685411299E-06
               0.1352460571822227E-07 -0.2969583000777938E-07  0.6374265424413608E+00  0.4391919124687460E-01  0.0000000000000000E+00  0.5592919170820314E-14
              -0.4079601917518316E-05  0.1052556482124766E-05 -0.1351773220872402E+02  0.6374258159142916E+00  0.0000000000000000E+00 -0.4735132312703365E-13
               0.1964238541540386E-06 -0.3720806396318763E-07 -0.8403098936929155E-14 -0.1661103086541925E-15  0.1000000000000000E+01 -0.2358669003370585E+01
               0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  ]

    @test GTPSA.jacobian(b0.coords.v) ≈ M_ESR

    p0c = 10e6
    # E to R_ref
    R_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + BeamTracking.massof(Species("electron"))^2))

    # Thin straight pure dipole:
    ele = LineElement(L=0.0, Kn0L=0.1, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin straight dipole:
    ele = LineElement(L=0.0, Kn0L=0.1, Kn5L=1000.0, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin pure quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, Kn5L=100.0, tilt5=-0.1*pi, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin pure multipole:
    ele = LineElement(L=0.0, Kn3L=10.0, tilt3=0.5*pi, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin multipole:
    ele = LineElement(L=0.0, Kn2L=1.0, tilt2=0.3*pi, Kn6L=100.0, tilt6=0.15*pi, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Drift: 
    ele = LineElement(L=1.0, tracking_method=SplitIntegration())   
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected = read_map("bmad_maps/drift.jl")
    q_expected = Quaternion(TPS64{D10}(1), TPS64{D10}[0, 0, 0])
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 0.0)

    # Curved drift: 
    ele = LineElement(L=2.0, g_ref=0.1, tracking_method=SplitIntegration())   
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/bend_no_field.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # Rotated and curved drift: 
    ele = LineElement(L=2.0, g_ref=0.1, tilt_ref=-pi/3, tracking_method=SplitIntegration())   
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/rotated_bend_no_field.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # Pure bend:
    ele = LineElement(L=2.0, g=0.1, tracking_method=SplitIntegration(order=6, num_steps=10))   
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    @test b0.coords.v ≈ [58.61782947 31.13531470 105.79375452 40.00000000 41.75205992 60.00000000]/1e3
    @test b0.coords.q ≈ [0.99999555887473 0.00000197685011 0.00297918168991 0.00008187412527]

    # Pure solenoid:
    ele = LineElement(L=1.0, Ksol=2.0, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Solenoid with quadrupole:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sol_quad.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # SK multiple steps:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=4, num_steps=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Step size:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=4, ds_step=1.0))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Straight pure dipole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, tilt0=pi/3, tracking_method=DriftKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_pure_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Straight pure dipole (BK):
    ele = LineElement(L=2.0, Kn0=0.1, tracking_method=BendKick(order=6, num_steps=10))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_pure_dipole_bk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-6)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-6)

    # Straight dipole with quadrupole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.03, tracking_method=DriftKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Straight dipole with quadrupole (BK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.1, tracking_method=BendKick(order=6, num_steps=10))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_dipole_bk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-6)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-6)

    # Straight dipole with quadrupole (MK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.1, tracking_method=SplitIntegration(order=6, num_steps=10))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_dipole_bk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 4e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 4e-7)

    # Pure quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quadrupole_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Skew quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=pi/4, tracking_method=MatrixKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-9)

    # Skew quadrupole another way (MK):
    ele = LineElement(L=2.0, Ks1=-0.1, tracking_method=MatrixKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-9)

    # Pure quadrupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=0.1*pi, tracking_method=DriftKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quadrupole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Quadrupole with octupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quad_oct_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # MK multiple steps:
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(num_steps=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/mk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Quadrupole with octupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=DriftKick())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quad_oct_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 4e-8)

    # Pure sextupole:
    ele = LineElement(L=2.0, Kn2=10.0, tilt2=0.2*pi, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sextupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 7e-8)

    # DK multiple steps:
    ele = LineElement(L=2.0, Kn2=10.0, tracking_method=SplitIntegration(order=4, num_steps=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/dk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Sextupole with decapole:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sex_dec.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-8)

    # Order 4:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=4))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_four.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Order 6:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=6))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_six.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-8)

    # Order 8:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=SplitIntegration(order=8))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_eight.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 3e-7)

    # Patch:
    ele = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=SplitIntegration())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # Particle lost in dipole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_dipole = LineElement(L=1.0, Kn0=1e-8, Kn1=1e-8, tracking_method=BendKick())
    track!(b0, Beamline([ele_dipole], R_ref=R_ref))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Particle lost in quadrupole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_quad = LineElement(L=1.0, Kn1=1e-8, tracking_method=MatrixKick())
    track!(b0, Beamline([ele_quad], R_ref=R_ref))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Particle lost in patch (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=SplitIntegration())
    track!(b0, Beamline([ele_patch], R_ref=R_ref))
    @test b0.coords.state[1] == State.Lost
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Errors:
    @test_throws ErrorException MatrixKick(ds_step = 0.1, num_steps = 2)
    @test_throws ErrorException BendKick(order = 2, num_steps = -2)
    @test_throws ErrorException DriftKick(ds_step = -0.1)
    @test_throws ErrorException SolenoidKick(num_steps = -2)
    @test_throws ErrorException SplitIntegration(order = 5)
  end

end