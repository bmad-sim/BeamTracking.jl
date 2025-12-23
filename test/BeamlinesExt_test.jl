@testset "Beamlines" begin

  include("BeamlinesExt/beamlines_utils_test.jl")
  include("BeamlinesExt/beamlines_aperture_test.jl")
  include("BeamlinesExt/beamlines_alignment_test.jl")
  include("BeamlinesExt/beamlines_stochastic_test.jl")

  #------------------------------------------------------------------------------------------------

  include("lattices/esr.jl")

  @testset "Exact" begin
    p0c = 10e6
    # E to R_ref
    R_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + massof(Species("electron"))^2))
    
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
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost (abs(px) > pt):
    b0 = Bunch([0.1 0.2 0.3 0.4 0.5 -0.6], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in bend (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_bend = LineElement(L=1.0, g_ref=0.01, Kn0=0.01, tracking_method=Exact())
    track!(b0, Beamline([ele_bend], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in drift (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_drift = LineElement(L=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_drift], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v

    # Particle lost in solenoid (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    ele_sol = LineElement(L=1.0, Ksol=1.0, tracking_method=Exact())
    track!(b0, Beamline([ele_sol], R_ref=R_ref, species_ref=Species("electron")))
    @test b0.coords.state[1] == STATE_LOST
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
  
  @testset "Yoshida" begin
    b0 = Bunch(collect(transpose(@vars(D1))), R_ref=ring.R_ref)
    foreach(t->t.tracking_method=Yoshida(), ring.line)
    track!(b0, ring)
    M_ESR = [0.8763088913632153E+00  0.2842332738570844E+00 -0.9233408564026070E-06 -0.1104742929395010E-06  0.1231581327803036E-08 -0.8939291467979220E-07 
            -0.8165279324836996E+00  0.8763069736287663E+00 -0.1898521122770908E-05 -0.1113630178379456E-06 -0.2820752660094096E-08 -0.1395543922780640E-06
             0.1352460499137301E-07 -0.2969583027665362E-07  0.6374265424413070E+00  0.4391919124687560E-01 -0.3331954964206950E-16  0.2618845432086097E-14 
            -0.4079601894069816E-05  0.1052556484633936E-05 -0.1351773220872471E+02  0.6374258159142887E+00 -0.4040958440539453E-14 -0.3122719464052644E-13
             0.1858496620330223E-06 -0.3179065075865608E-07  0.1927488251425578E-13 -0.2937455989329897E-14  0.9343600246296299E+00 -0.2307665523810159E+01  
             0.6685543985517410E-08 -0.3425164613573595E-08  0.2907237074330440E-14  0.1826161170763475E-15  0.4150093714505364E-01  0.9677530013155135E+00]

    @test GTPSA.jacobian(b0.coords.v) ≈ M_ESR

    p0c = 10e6
    # E to R_ref
    R_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + massof(Species("electron"))^2))

    # Thin straight pure dipole:
    ele = LineElement(L=0.0, Kn0L=0.1, tracking_method=Yoshida())
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
    ele = LineElement(L=0.0, Kn0L=0.1, Kn5L=1000.0, tracking_method=Yoshida())
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
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, tracking_method=Yoshida())
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
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, Kn5L=100.0, tilt5=-0.1*pi, tracking_method=Yoshida())
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
    ele = LineElement(L=0.0, Kn3L=10.0, tilt3=0.5*pi, tracking_method=Yoshida())
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
    ele = LineElement(L=0.0, Kn2L=1.0, tilt2=0.3*pi, Kn6L=100.0, tilt6=0.15*pi, tracking_method=Yoshida())
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
    ele = LineElement(L=1.0, tracking_method=Yoshida())   
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

    # Drift another way: 
    ele = LineElement(L=1.0, Kn0=0.0, Kn1=0.0, tracking_method=MatrixKick())   
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

    # Drift yet another way: 
    ele = LineElement(L=1.0, Kn2=0.0, Kn1=0.0, tracking_method=MatrixKick())   
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

    # Drift yet yet another way: 
    ele = LineElement(L=1.0, Kn1=0.0, tracking_method=MatrixKick())   
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
    ele = LineElement(L=2.0, g_ref=0.1, tracking_method=Yoshida())   
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
    ele = LineElement(L=2.0, g_ref=0.1, tilt_ref=-pi/3, tracking_method=Yoshida())   
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
    ele = LineElement(L=2.0, g=0.1, tracking_method=Yoshida(order=6, num_steps=10))   
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    @test b0.coords.v ≈ [58.61782947 31.13531470 105.79375452 40.00000000 41.75205992 60.00000000]/1e3
    @test (b0.coords.q ≈ [0.99999555887473 0.00000197685011 0.00297918168991 0.00008187412527]
           || b0.coords.q ≈ -[0.99999555887473 0.00000197685011 0.00297918168991 0.00008187412527])

    # Pure solenoid:
    ele = LineElement(L=1.0, Ksol=2.0, tracking_method=Yoshida(order=2))
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
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=2))
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
    ele = LineElement(L=2.0, Kn0=0.1, tilt0=pi/3, tracking_method=DriftKick(order=2))
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
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.03, tracking_method=DriftKick(order=2))
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
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.1, tracking_method=Yoshida(order=6, num_steps=10))
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
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(order=2))
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
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=pi/4, tracking_method=MatrixKick(order=2))
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
    ele = LineElement(L=2.0, Ks1=-0.1, tracking_method=MatrixKick(order=2))
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
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=0.1*pi, tracking_method=DriftKick(order=2))
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
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(order=2))
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
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(order=2, num_steps=2))
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
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=DriftKick(order=2))
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
    ele = LineElement(L=2.0, Kn2=10.0, tilt2=0.2*pi, tracking_method=Yoshida(order=2))
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
    ele = LineElement(L=2.0, Kn2=10.0, tracking_method=Yoshida(order=4, num_steps=2))
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
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=2))
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
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=4))
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
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=6))
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
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=8))
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
    ele = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # RF Cavity (PTC):
    ele = LineElement(L=4.01667, voltage=3.3210942126011E6, rf_frequency=5.9114268014977E8, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/pure_rf.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # Harmon:
    ele_drift = LineElement(L=1.04812778909)
    ele = LineElement(L=4.01667, voltage=3.3210942126011E6, harmon=10, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele_drift, ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl.line[2])
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/pure_rf.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # With solenoid (RK4):
    ele = LineElement(L=4.01667, voltage=3321.0942126011, rf_frequency=591142.68014977, Ksol=0.6, tracking_method=Yoshida(order=6, num_steps=2))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.14844043934151196 0.01019264367221814 -0.00269118775927293 -0.00153213180245353 0.046619484431963405 0.0600001992395169]
    q_expected = [0.4182448759327037 0.00112479589051007 -0.00026561236717118287 -0.9083335775145173]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected || b0.coords.q ≈ -q_expected

    # With sextupole:
    ele = LineElement(L=4.01667, voltage=3321.0942126011, rf_frequency=591142.68014977, phi0=0.1, Kn2=1.3, tracking_method=Yoshida(order=6, num_steps=20))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.11924424064353095 0.05139343591364565 0.22257164275272448 0.08178778488624713 0.04410474005001025 0.05996700254315476]
    q_expected = [0.9996797773832821 -0.020232635960537492 0.015198115217757237 -2.0659959944872022e-5]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected || b0.coords.q ≈ -q_expected

    # With solenoid and quadrupole:
    ele = LineElement(L=4.01667, voltage=3321.0942126011, rf_frequency=591142.68014977, Ksol=-0.3, Kn1=0.15, tracking_method=Yoshida(order=6, num_steps=20))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, R_ref=R_ref, species=Species("electron"))
    bl = Beamline([ele], R_ref=R_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.06425313141616465 -0.016053284947290265 0.2828816007928012 0.11181265178119304 0.04054713815327194 0.06000019363164233]
    q_expected = [0.8408669037648832 -0.03581645006820144 0.0065436536214155076 0.5400159374080092]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected || b0.coords.q ≈ -q_expected

    # Bend with deterministic radiation:
    ele = LineElement(L=2.0, g=0.1, tracking_method=BendKick(order=6, radiation_damping_on=true))
    v = zeros(6)
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.0001092610284174973 -0.00016340536870756257 0.0 0.0 5.461341618518875e-6 -0.0016395105602759578]
    @test b0.coords.v ≈ v_expected

    # Quadrupole with deterministic radiation:
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(order=6, num_steps = 2, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.043612346633677884 0.014463450553206672 0.11624841932128559 0.05417994410649243 0.047841511069009836 0.05998800025940442]
    @test b0.coords.v ≈ v_expected

    # Sextupole with deterministic radiation:
    ele = LineElement(L=0.8, Kn2=1.3, tracking_method=DriftKick(order=6, num_steps = 5, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.025386841257333655 0.02092955214032272 0.06046374439640514 0.04086917032430656 0.04927228624504675 0.059999784343950334]
    @test b0.coords.v ≈ v_expected

    # Solenoid with deterministic radiation:
    ele = LineElement(L=1.5, Ksol=0.3, tracking_method=SolenoidKick(order=6, num_steps = 2, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.055079725897422445 0.026845723521363107 0.07564278174163326 0.03323733870836322 0.04860800161184006 0.05997689477653424]
    @test b0.coords.v ≈ v_expected

    # Cavity-solenoid with deterministic radiation:
    ele = LineElement(L=0.5, Ksol=0.3, rf_frequency=1e8, voltage=-0.25e6, tracking_method=Yoshida(order=6, num_steps = 5, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.022813869607319508 0.02259460588210925 0.04729868720097375 0.038077653342022774 0.04953599990242161 0.05999085207779174]
    @test b0.coords.v ≈ v_expected

    # Map:
    function transport_map(v, q)
      v_out = (sin(v[1]), 2*v[2], exp(v[3]), 1-v[4], v[6], v[5])
      if !isnothing(q)
        q_out = (q[2], q[1], q[4], q[3])
      else
        q_out = nothing
      end
      return v_out, q_out
    end

    ele = LineElement(transport_map=transport_map)
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [sin(0.01) 2*0.02 exp(0.03) 1-0.04 0.06 0.05]
    @test b0.coords.v ≈ v_expected

    ele = LineElement(transport_map=transport_map)
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1/sqrt(2) 0.0 1/sqrt(2) 0.0]
    b0 = Bunch(v, q, R_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], R_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [sin(0.01) 2*0.02 exp(0.03) 1-0.04 0.06 0.05]
    q_expected = [0.0 1/sqrt(2) 0.0 1/sqrt(2)]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected || b0.coords.q ≈ -q_expected

    # Particle lost in dipole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_dipole = LineElement(L=1.0, Kn0=1e-8, Kn1=1e-8, tracking_method=BendKick())
    track!(b0, Beamline([ele_dipole], R_ref=R_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Particle lost in quadrupole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_quad = LineElement(L=1.0, Kn1=1e-8, tracking_method=MatrixKick())
    track!(b0, Beamline([ele_quad], R_ref=R_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Particle lost in patch (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], R_ref=R_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Yoshida())
    track!(b0, Beamline([ele_patch], R_ref=R_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q || q_init == -b0.coords.q

    # Errors:
    @test_throws ErrorException MatrixKick(ds_step = 0.1, num_steps = 2)
    @test_throws ErrorException BendKick(order = 2, num_steps = -2)
    @test_throws ErrorException DriftKick(ds_step = -0.1)
    @test_throws ErrorException SolenoidKick(num_steps = -2)
    @test_throws ErrorException Yoshida(order = 5)
  end  
end