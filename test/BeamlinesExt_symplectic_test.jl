@testset "Beamlines" begin

  include("lattices/esr.jl")

  @testset "Yoshida" begin
    b0 = Bunch(collect(transpose(@vars(D1))), p_over_q_ref=ring.p_over_q_ref)
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
    # E to p_over_q_ref
    p_over_q_ref = BeamTracking.E_to_R(Species("electron"), sqrt(p0c^2 + massof(Species("electron"))^2))

    # Thin straight pure dipole:
    ele = LineElement(L=0.0, Kn0L=0.1, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin straight dipole:
    ele = LineElement(L=0.0, Kn0L=0.1, Kn5L=1000.0, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_dipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin pure quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin quadrupole:
    ele = LineElement(L=0.0, Kn1L=1.0, tilt1=pi, Kn5L=100.0, tilt5=-0.1*pi, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_quadrupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin pure multipole:
    ele = LineElement(L=0.0, Kn3L=10.0, tilt3=0.5*pi, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_pure_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Thin multipole:
    ele = LineElement(L=0.0, Kn2L=1.0, tilt2=0.3*pi, Kn6L=100.0, tilt6=0.15*pi, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/thin_multipole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 1e-14)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Drift: 
    ele = LineElement(L=1.0, tracking_method=Yoshida())   
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
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
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/bend_no_field.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # Rotated and curved drift: 
    ele = LineElement(L=2.0, g_ref=0.1, tilt_ref=-pi/3, tracking_method=Yoshida())   
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/rotated_bend_no_field.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # Pure bend:
    ele = LineElement(L=2.0, g=0.1, tracking_method=Yoshida(order=6, n_steps=10, fringe_at=Fringe.NoEnd))   
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    @test b0.coords.v ≈ [58.61782947 31.13531470 105.79375452 40.00000000 41.75205992 60.00000000]/1e3
    @test b0.coords.q ≈ [0.99999555887473 0.00000197685011 0.00297918168991 0.00008187412527]

    # Pure solenoid:
    ele = LineElement(L=1.0, Ksol=2.0, tracking_method=Yoshida(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected, _ = read_spin_orbit_map("bmad_maps/solenoid.jl")
    q_expected = [0.0                   0.0                   0.0                  0.0                  0.0 0.8430735212133056; 
                 -0.018145294992663086  0.011650959928739194  0.011650959928739194 0.018145294992663097 0.0 0.0; 
                 -0.011650959928739215 -0.018145294992663066 -0.018145294992663066 0.011650959928739194 0.0 0.0; 
                  0.0                   0.0                   0.0                  0.0                  0.0 0.5399515598487606]
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test scalar.(b0.coords.q) ≈ [0.5393261291271401 0.0 0.0 -0.8420969816124171]
    @test GTPSA.jacobian(b0.coords.q) ≈ q_expected

    # Solenoid with quadrupole:
    ele = LineElement(L=2.0, Ksol=0.1, Kn1=0.1, tracking_method=SolenoidKick(order=6, n_steps=8, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06] .+ collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.8084528550552892      1.7618358451908693    0.10016574979271514    0.1799542070093867  0.0 -0.045295756768258456; 
                 -0.1914764787788115      0.8086450231845644   -0.00043899255470597357 0.08804905786594175 0.0  0.0005199878910795101; 
                 -0.08837437582642556    -0.17667829248457204   1.1858656497844702     2.00204600773925    0.0 -0.07705147167235063; 
                  0.00043023645106843864 -0.09983080426545739   0.20723515702913622    1.1856282544251113  0.0 -0.005074322385922071; 
                  0.007771095889150356   -0.030748733815372683 -0.009982509694960352  -0.08527719770242106 1.0  0.008496403081924352; 
                  0.0                     0.0                   0.0                    0.0                 0.0  1.0]
    q_expected = [-0.0001968320381634533 -0.0003354522139518879 -0.0007323964953026688  -0.0010197270511281822 0.0 0.008517335358001279; 
                   0.003945765503469646   0.010191016877089102  -0.102261685137679      -0.09380414899291793   0.0 0.009442612029680206; 
                  -0.09027034467726167   -0.08810488645733103   -0.00479577073391418    -0.006822247275792119  0.0 0.005170996125814023; 
                   0.0005861140247321611 -0.0013797021734982626 -0.00046618009310879727 -0.004017421548244223  0.0 0.08899043520023708]
    @test scalar.(b0.coords.v) ≈ [0.05342292167847618 0.017767436296301228 0.11106668172142604 0.05163945550900599 0.04816617224636819 0.06]
    @test scalar.(b0.coords.q) ≈ [0.9955073952270296 -0.0065554530183620785 -0.003070190160659054 -0.09440670535724088]
    @test GTPSA.jacobian(b0.coords.v) ≈ v_expected
    @test GTPSA.jacobian(b0.coords.q) ≈ q_expected

    # Straight pure dipole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, tilt0=pi/3, tracking_method=DriftKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_pure_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Straight pure dipole (BK):
    ele = LineElement(L=2.0, Kn0=0.1, tracking_method=BendKick(order=2, n_steps=1, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_pure_dipole_bk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-6)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-6)

    # Straight dipole with quadrupole (DK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.03, tracking_method=DriftKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_dipole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Straight dipole with quadrupole (BK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.1, tracking_method=BendKick(order=6, n_steps=10, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06] .+ collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.8134564588420554     1.7929903599772967    -0.0009771147239118648 -0.007170466517246891  0.0  0.14076980478413328; 
                 -0.18750993419337825    0.8160197842435812    4.353810553546985e-5    0.0003895479880845047 0.0 -0.008330646843440145; 
                  0.0009307219466424009 -0.006066641914690669  1.1966220762978186      2.022908831874305     0.0 -0.08753711731617767; 
                  4.207781984425512e-5  -0.0003478113023711722 0.21290312674959588     1.195601946695064     0.0 -0.008018012610111421; 
                 -0.019615339181106393   0.1297895191980294   -0.009044406233040323   -0.08844473737824309   1.0  0.023971876629362485; 
                  0.0 0.0 0.0 0.0 0.0 1.0]
    q_expected = [-0.008604423440626487  -0.007827018117037557  -0.0007491418589629741 -0.0010109567445952526 0.0  0.007860003014201084; 
                  -0.00010027511069209182 0.0007838650281191009 -0.1032835902104188    -0.09551631148916687   0.0  0.010323329616161184; 
                  -0.09149168223931535   -0.08329012934091672   -0.0003547892058559101 -0.003700238847439373  0.0  0.08282160471534457; 
                   0.000222787405508111   0.00012056511240449336 0.0005709484656144937  0.005181887546895218  0.0 -0.00048747250073165424]
    @test scalar.(b0.coords.v) ≈ [-0.14043738611320963 -0.17314263313440387 0.1166396952761015 0.054196794002912455 0.03990991202603635 0.06]
    @test scalar.(b0.coords.q) ≈ [0.9955838402416232 -0.006898381615103346 -0.09362253681668872 0.0002235639518493756]
    @test GTPSA.jacobian(b0.coords.v) ≈ v_expected
    @test GTPSA.jacobian(b0.coords.q) ≈ q_expected

    # Straight dipole with quadrupole (MK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.1, tracking_method=Yoshida(order=6, n_steps=10, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/straight_dipole_bk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 4e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 4e-7)

    # Pure quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quadrupole_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Skew quadrupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=pi/4, tracking_method=MatrixKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-9)

    # Skew quadrupole another way (MK):
    ele = LineElement(L=2.0, Ks1=-0.1, tracking_method=MatrixKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/skew_quad_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-9)

    # Pure quadrupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, tilt1=0.1*pi, tracking_method=DriftKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quadrupole_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-9)

    # Quadrupole with octupole (MK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quad_oct_mk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # MK multiple steps:
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=MatrixKick(order=2, n_steps=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/mk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Quadrupole with dipole and sextupole (MK):
    ele = LineElement(L=2.0, Kn0=0.1, Kn1=0.2, Kn2=0.3, tracking_method=MatrixKick(order=6, n_steps=10, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.13784912725589 -0.16544887326705 0.12734512769155 0.06742027838098 0.03975920093905 0.06]
    q_expected = [0.99586268249262 -0.01331711276472 -0.08988910144266 0.00034866609706]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected

    # Quadrupole with octupole (DK):
    ele = LineElement(L=2.0, Kn1=0.1, Kn3=100.0, tracking_method=DriftKick(order=2, fringe_at=Fringe.NoEnd))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/quad_oct_dk.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 4e-8)

    # Pure sextupole:
    ele = LineElement(L=2.0, Kn2=10.0, tilt2=0.2*pi, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sextupole.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 7e-8)

    # DK multiple steps:
    ele = LineElement(L=2.0, Kn2=10.0, tracking_method=Yoshida(order=4, n_steps=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/dk_multistep.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Sextupole with decapole:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/sex_dec.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-8)

    # Order 4:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=4))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_four.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 2e-7)

    # Order 6:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=6))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_six.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 6e-8)

    # Order 8:
    ele = LineElement(L=2.0, Kn2=10.0, Kn4=100.0, tracking_method=Yoshida(order=8))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/order_eight.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-9)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 3e-7)

    # Patch:
    ele = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Yoshida())
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/patch.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 5e-10)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-14)

    # RF Cavity (PTC):
    ele = LineElement(L=4.01667, voltage=3.3210942126011E6, zero_phase=PhaseRef.AboveTransition, rf_frequency=5.9114268014977E8, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/pure_rf.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # Harmon:
    ele_drift = LineElement(L=1.04812778909)
    ele = LineElement(L=4.01667, voltage=3.3210942126011E6, harmon=10, zero_phase=PhaseRef.AboveTransition, tracking_method=Yoshida(order=2))
    v = collect(transpose(@vars(D10)))
    q = TPS64{D10}[1 0 0 0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele_drift, ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl.line[2])
    q_z = Quaternion(b0.coords.q[1], b0.coords.q[2:4])
    v_expected, q_expected = read_spin_orbit_map("bmad_maps/pure_rf.jl")
    @test coeffs_approx_equal(v_expected, b0.coords.v, 2e-7)
    @test quaternion_coeffs_approx_equal(q_expected, q_z, 1e-7)

    # With solenoid (RK4):
    ele = LineElement(L=4.01667, voltage=3321.0942126011,  zero_phase=PhaseRef.AboveTransition, rf_frequency=591142.68014977, Ksol=0.6, tracking_method=Yoshida(order=6, n_steps=2, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.14844043934151196 0.01019264367221814 -0.00269118775927293 -0.00153213180245353 0.046619484431963405 0.0600001992395169]
    q_expected = [0.4182448759327037 0.00112479589051007 -0.00026561236717118287 -0.9083335775145173]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected

    # With sextupole:
    ele = LineElement(L=4.01667, voltage=3321.0942126011, zero_phase=PhaseRef.AboveTransition,  rf_frequency=591142.68014977, phi0=0.1, Kn2=1.3, tracking_method=Yoshida(order=6, n_steps=20))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.11924424064353095 0.05139343591364565 0.22257164275272448 0.08178778488624713 0.04410474005001025 0.05996700254315476]
    q_expected = [0.9996797773832821 -0.020232635960537492 0.015198115217757237 -2.0659959944872022e-5]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected

    # With solenoid and quadrupole:
    ele = LineElement(L=4.01667, voltage=3321.0942126011,  zero_phase=PhaseRef.AboveTransition, rf_frequency=591142.68014977, Ksol=-0.3, Kn1=0.15, tracking_method=Yoshida(order=6, n_steps=20, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.06425313141616465 -0.016053284947290265 0.2828816007928012 0.11181265178119304 0.04054713815327194 0.06000019363164233]
    q_expected = [0.8408669037648832 -0.03581645006820144 0.0065436536214155076 0.5400159374080092]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected

    # Bend with deterministic radiation:
    ele = LineElement(L=2.0, g=0.1, tracking_method=BendKick(order=6, radiation_damping_on=true))
    v = zeros(6)
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.0001092610284174973 -0.00016340536870756257 0.0 0.0 5.461341618518875e-6 -0.0016395105602759578]
    @test b0.coords.v ≈ v_expected

    # Quadrupole with deterministic radiation:
    ele = LineElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(order=6, n_steps=2, radiation_damping_on=true, fringe_at=Fringe.NoEnd))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.043612346633677884 0.014463450553206672 0.11624841932128559 0.05417994410649243 0.047841511069009836 0.05998800025940442]
    @test b0.coords.v ≈ v_expected

    # Sextupole with deterministic radiation:
    ele = LineElement(L=0.8, Kn2=1.3, tracking_method=DriftKick(order=6, n_steps=5, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.025386841257333655 0.02092955214032272 0.06046374439640514 0.04086917032430656 0.04927228624504675 0.059999784343950334]
    @test b0.coords.v ≈ v_expected

    # Solenoid with deterministic radiation:
    ele = LineElement(L=1.5, Ksol=0.3, tracking_method=SolenoidKick(order=6, n_steps=2, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.055079725897422445 0.026845723521363107 0.07564278174163326 0.03323733870836322 0.04860800161184006 0.05997689477653424]
    @test b0.coords.v ≈ v_expected

    # Cavity-solenoid with deterministic radiation:
    ele = LineElement(L=0.5, Ksol=0.3, rf_frequency=1e8,  zero_phase=PhaseRef.AboveTransition, voltage=-0.25e6, tracking_method=Yoshida(order=6, n_steps = 5, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
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
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [sin(0.01) 2*0.02 exp(0.03) 1-0.04 0.06 0.05]
    @test b0.coords.v ≈ v_expected

    ele = LineElement(transport_map=transport_map)
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1/sqrt(2) 0.0 1/sqrt(2) 0.0]
    b0 = Bunch(v, q, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [sin(0.01) 2*0.02 exp(0.03) 1-0.04 0.06 0.05]
    q_expected = [0.0 1/sqrt(2) 0.0 1/sqrt(2)]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected

    # Implicit:
    function crazy(x, y, s, t, p)
      a0, c0 = p
      potential = (a0*C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), a0*cosh(x)*sin(y)*cos(s)*sinh(t), a0*sinh(x)*cosh(y)*sin(s)*cos(t), a0*c0*cos(x)*sinh(y)*cosh(s)*sin(t))
      jac =  (a0*C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -a0*C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), a0*C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), a0*C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
          a0*cos(s)*sin(y)*sinh(t)*sinh(x),  a0*cos(s)*cos(y)*cosh(x)*sinh(t), -a0*cosh(x)*sin(s)*sin(y)*sinh(t),  a0*cos(s)*cosh(t)*cosh(x)*sin(y),
          a0*cos(t)*cosh(x)*cosh(y)*sin(s),  a0*cos(t)*sin(s)*sinh(x)*sinh(y),  a0*cos(s)*cos(t)*cosh(y)*sinh(x), -a0*cosh(y)*sin(s)*sin(t)*sinh(x),
        -a0*c0*cosh(s)*sin(t)*sin(x)*sinh(y), a0*c0*cos(x)*cosh(s)*cosh(y)*sin(t), a0*c0*cos(x)*sin(t)*sinh(s)*sinh(y), a0*c0*cos(t)*cos(x)*cosh(s)*sinh(y))
      return potential, jac
    end
    ele = LineElement(four_potential=crazy, four_potential_params=(1.0, 1.0), four_potential_normalized=true, L=1.2, tracking_method=Yoshida(order=8, n_steps=60, radiation_damping_on=true))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.24540112246103782 -0.7846013184552288 0.12024259722147093 0.031368656927501536 -0.004848169943849862 0.05999999971081322]
    q_expected = [0.9149630270850545 -0.0023107478986964395 -0.31231391206369696 -0.2555334417275888]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected 

    function dipole(x, y, s, t, p)
      g, Kn0 = p
       potential = (0.0, 0.0, 0.0, -Kn0*(x + g*x^2/2))
        jac = (0.0,            0.0, 0.0, 0.0,
               0.0,            0.0, 0.0, 0.0,
               0.0,            0.0, 0.0, 0.0,
               -Kn0*(1 + g*x), 0.0, 0.0, 0.0)
      return potential, jac
    end
    ele = LineElement(L=2.0, g_ref=0.1, four_potential=dipole, four_potential_params=(0.1, 0.1), four_potential_normalized=true, tracking_method=Yoshida(order=6, n_steps=10, radiation_damping_on=true))
    v = zeros(6)
    b0 = Bunch(v, p_over_q_ref=-18e9/C_LIGHT, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=-18e9/C_LIGHT, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.0001092610284174973 -0.00016340536870756257 0.0 0.0 5.461341618518875e-6 -0.0016395105602759578]
    @test b0.coords.v ≈ v_expected

    # Solenoid hard-edge fringe:
    ele = LineElement(Ksol=2.0, L=1.0, tracking_method=SolenoidKick(order=2, fringe_at=Fringe.BothEnds))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.05344208147045112 0.0014068714384027799 0.011406871438402784 -0.0034420814704511204 0.0486268409937747 0.06]
    q_expected = [0.5850869992845447 0.005047034551437897 -0.020079160216966062 -0.8107062094466929]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected 

    # Dipole hard-edge fringe:
    ele = LineElement(Kn0=0.1, L=1.0, tracking_method=BendKick(order=2, fringe_at=Fringe.BothEnds))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [-0.0185407314436441 -0.07999999999999999 0.06773781043103988 0.0394311017804856 0.0486403579195499 0.06]
    q_expected = [0.9988283540819322 0.00022881568956878987 -0.04835817826999246 0.0018312071881943616]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected 

    # Quadrupole hard-edge fringe:
    ele = LineElement(Kn1=0.36, L=1.0, tracking_method=MatrixKick(order=6, n_steps=10, fringe_at=Fringe.BothEnds))
    v = [0.01 0.02 0.03 0.04 0.05 0.06]
    q = [1.0 0.0 0.0 0.0]
    b0 = Bunch(v, q, p_over_q_ref=p_over_q_ref, species=Species("electron"))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    track!(b0, bl)
    v_expected = [0.026172709210547637 0.013283654892913688 0.0752222208993626 0.0583923740672597 0.048976075148431955 0.06]
    q_expected = [0.9999550310129977 -0.008905896376194768 -0.0032459810037083338 0.00029080723839939075]
    @test b0.coords.v ≈ v_expected
    @test b0.coords.q ≈ q_expected 

    # Particle lost in dipole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_dipole = LineElement(L=1.0, Kn0=1e-8, Kn1=1e-8, tracking_method=BendKick(fringe_at=Fringe.NoEnd))
    track!(b0, Beamline([ele_dipole], p_over_q_ref=p_over_q_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q

    # Particle lost in quadrupole (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_quad = LineElement(L=1.0, Kn1=1e-8, tracking_method=MatrixKick(fringe_at=Fringe.NoEnd))
    track!(b0, Beamline([ele_quad], p_over_q_ref=p_over_q_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q

    # Particle lost in patch (momentum is too small):
    b0 = Bunch([0.4 0.4 0.4 0.4 0.4 -0.5], [1.0 0.0 0.0 0.0], p_over_q_ref=p_over_q_ref, species=Species("electron"))
    v_init = copy(b0.coords.v)
    q_init = copy(b0.coords.q)
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Yoshida())
    track!(b0, Beamline([ele_patch], p_over_q_ref=p_over_q_ref))
    @test b0.coords.state[1] == STATE_LOST
    @test v_init == b0.coords.v
    @test q_init == b0.coords.q

    # Errors:
    @test_throws ErrorException MatrixKick(ds_step = 0.1, n_steps = 2)
    @test_throws ErrorException BendKick(order = 2, n_steps = -2)
    @test_throws ErrorException DriftKick(ds_step = -0.1)
    @test_throws ErrorException SolenoidKick(n_steps = -2)
    @test_throws ErrorException Yoshida(order = 5)
  end  
end
