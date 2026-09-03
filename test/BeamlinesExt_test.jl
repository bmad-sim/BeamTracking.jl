@testset "Beamlines" begin

  include("BeamlinesExt/beamlines_sagan_cavity_test.jl")
  include("BeamlinesExt/beamlines_utils_test.jl")
  include("BeamlinesExt/beamlines_aperture_test.jl")
  include("BeamlinesExt/beamlines_alignment_test.jl")
  include("BeamlinesExt/beamlines_stochastic_test.jl")
  include("BeamlinesExt/beamlines_ibs_test.jl")
  include("BeamlinesExt/beamlines_context_test.jl")

  #------------------------------------------------------------------------------------------------


  @testset "Exact" begin
    p0c = 10e6
    # E to p_over_q_ref
    p_over_q_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + massof(Species("electron"))^2))
    
    # Patch:
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_patch], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Drift: 
    ele_drift = LineElement(L=1.0, tracking_method=Exact())   
    b0 = Bunch(collect(transpose(@vars(D10))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_drift], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/drift.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Thick solenoid:
    ele_sol = LineElement(L=1.0, Ksol=2.0, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D10))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_sol], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = read_map("bmad_maps/solenoid.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)

    # Thick pure bend:
    ele_bend = LineElement(L=0.5, g_ref=1, Kn0 = 1.001, tracking_method=Exact())
    b0 = Bunch(collect(transpose(@vars(D1))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_thick_dipole], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps2)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps3)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps4)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps5)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps6)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps7)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps8)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(collect(transpose(@vars(D1) + ps9)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)

    @test GTPSA.jacobian(b0.coords.v) ≈ exact_bend_15

    # Particle lost (does not intersect exit face):
    b0 = Bunch(collect(transpose(@vars(D1))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_kick = LineElement(L=1.0, Kn0L=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_kick], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost (abs(px) > pt):
    b0 = Bunch([0.1 0.2 0.3 0.4 0.5 -0.6], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in bend (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in drift (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_drift = LineElement(L=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_drift], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in solenoid (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_sol = LineElement(L=1.0, Ksol=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_sol], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Errors:
    b0 = Bunch(collect(transpose(@vars(D1))), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_patch_bend = LineElement(L=1.0, g_ref=0.01, dy=3.0, dz_rot=0.3, tracking_method=Exact())
    ele_patch_sol = LineElement(L=1.0, Ksol=1.0, dt=1.0, tracking_method=Exact())
    ele_bend_quad = LineElement(L=1.0, g_ref=0.01, Kn1=1.0, tracking_method=Exact())
    @test_throws ErrorException track!(b0, Beamline([ele_patch_bend], p_over_q_ref=p_over_q_ref, species_ref=Species("electron")))
    @test_throws ErrorException track!(b0, Beamline([ele_patch_sol],  p_over_q_ref=p_over_q_ref))
    @test_throws ErrorException track!(b0, Beamline([ele_bend_quad],  p_over_q_ref=p_over_q_ref))
  end
end