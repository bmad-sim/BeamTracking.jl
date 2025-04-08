using Test,
      BeamTracking,
      JET,
      BenchmarkTools,
      GTPSA

BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
BenchmarkTools.DEFAULT_PARAMETERS.evals = 2
  
# Soon we generalize this to test_map...
function test_matrix(ele, n_work, M_expected; type_stable=true, no_allocs=true, tol=1e-14, 
                                              beta_gamma_ref=1.0, species=Species("electron"))
    
  GTPSA.desc_current = Descriptor(6, 1) # 6 variables, order 1
  bunch = Bunch(beta_gamma_ref=beta_gamma_ref, species=species, gtpsa_map=true) 
  work = BeamTracking.get_work(bunch, Val(n_work))

  track!(bunch, ele, work=work)
  p = Particle(bunch, 1)

  # 1) Correctness
  @test norm(GTPSA.jacobian(p.v) - M_expected) < tol 
  # 2) Type stability
  if type_stable
    @test_opt track!(bunch, ele, work=work)
  end
  # 3) No Allocations
  if no_allocs
    @test @ballocated(track!($bunch, $ele, work=$work)) == 0 
  end
end

@testset "BeamTracking.jl" begin
  # Linear tests -----------------------  
  beta_gamma_ref = 1.0
  gamma_ref = sr_gamma(beta_gamma_ref)
  species = Species("electron")
  brho_ref = brho(massof(species), beta_gamma_ref, chargeof(species)) 
  tol = 1e-14 

  # Linear drift test:
    L_d = 5.0
    drift = Linear.Drift(L=L_d)
    M_drift_expected = [1.0  L_d  0.0  0.0  0.0  0.0;
                        0.0  1.0  0.0  0.0  0.0  0.0;
                        0.0  0.0  1.0  L_d  0.0  0.0;
                        0.0  0.0  0.0  1.0  0.0  0.0;
                        0.0  0.0  0.0  0.0  1.0  L_d/gamma_ref^2;
                        0.0  0.0  0.0  0.0  0.0  1.0] 
                        
    test_matrix(drift, 0, M_drift_expected, beta_gamma_ref=beta_gamma_ref, species=species)
    
  # Linear Quadrupole test----------------------- 
    
    # Focusing:
    L_q = 1.2
    Kn1 = 0.36
    Bn1 = Kn1*brho_ref
    qf = Linear.Quadrupole(Bn1=Bn1,L=L_q)
    M_qf_x = [cos(sqrt(Kn1)*L_q)            sincu(sqrt(Kn1)*L_q)*L_q;  
              -sqrt(Kn1)*sin(sqrt(Kn1)*L_q) cos(sqrt(Kn1)*L_q)      ;]
    M_qf_y = [cosh(sqrt(Kn1)*L_q)           sinhcu(sqrt(Kn1)*L_q)*L_q; 
              sqrt(Kn1)*sinh(sqrt(Kn1)*L_q) cosh(sqrt(Kn1)*L_q)      ;]
    M_qf_z = [1.0 L_q/gamma_ref^2;
              0.0 1.0            ;] 
    M_qf_expected = zeros(6,6)
    M_qf_expected[1:2,1:2] = M_qf_x
    M_qf_expected[3:4,3:4] = M_qf_y
    M_qf_expected[5:6,5:6] = M_qf_z

    test_matrix(qf, 1, M_qf_expected, beta_gamma_ref=beta_gamma_ref, species=species)

    # Defocusing:
    qd = Linear.Quadrupole(Bn1=-Bn1,L=L_q)
    M_qd_expected = zeros(6,6)
    M_qd_expected[1:2,1:2] = M_qf_y
    M_qd_expected[3:4,3:4] = M_qf_x
    M_qd_expected[5:6,5:6] = M_qf_z 

    test_matrix(qd, 1, M_qd_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  #Linear Solenoid test----------------------- 
    L_s =1.0
    S = 1.57
    Bs = S * brho_ref
    phi = S * L_s / 2
    so = Linear.Solenoid(L = L_s, Bs = Bs)
    M_solenoid_expected = [cos(phi)^2  sin(2*phi)/S  sin(2*phi)/2  sin(phi)^2*2/S   0.0  0.0;
                          -sin(2*phi)*S/4  cos(phi)^2 -S/2*sin(phi)^2  sin(2*phi)/2  0.0 0.0;
                          -sin(2*phi)/2  -2/S*sin(phi)^2 cos(phi)^2 sin(2*phi)/S  0.0 0.0;
                          S/2*sin(phi)^2  -sin(2*phi)/2 -S/4*sin(2*phi) cos(phi)^2 0.0 0.0;
                            0.0 0.0 0.0 0.0 1 L_s/gamma_ref^2;
                            0.0 0.0 0.0 0.0 0.0 1]
    
    test_matrix(so, 3, M_solenoid_expected, beta_gamma_ref=beta_gamma_ref, species=species)


  #= 
    #Linear S-Bend Test----------------------- 
    L_sb = 1.0

    #design g = actual g 
    gtot = 0.25
    g = 0.25
    B0 = gtot * brho_ref
    e1 = 0.20
    e2 = 0.20
    gL = g*L_sb
    sb1 = Linear.SBend(L = L_sb, B0 = B0, g = g, e1 = e1, e2 = e2)
    M_sb_expected = [
      1.019063687076115200  0.98961583701809175   0.00000000000000000    0.00000000000000000    0.00000000000000000   0.12435031315742107
      0.038894687087008446  1.01906368707611520   0.00000000000000000    0.00000000000000000    0.00000000000000000   0.25370572335343677  
      0.00000000000000000   0.00000000000000000    0.94932249112283185   1.00000000000000000    0.00000000000000000   0.00000000000000000   
      0.00000000000000000   0.00000000000000000   -0.098786807848340791  0.94932249112283185    0.00000000000000000   0.00000000000000000
      -0.2537057233534368  -0.12435031315742107    0.00000000000000000    0.00000000000000000   1.00000000000000000   0.48961583701809164
      0.00000000000     0.00000000000     0.00000000000     0.00000000000     0.00000000000     1.00000000000
      ]

    test_matrix(sb1, 1, M_sb_expected, beta_gamma_ref=beta_gamma_ref, species=species)
    
    #design g unequal to actual g, no e 
    g = 0.249
    e1 = 0.00
    e2 = 0.00
    gL = sqrt(g*gtot)*L_sb
    sb2 = Linear.SBend(L = L_sb, B0 = B0, g = g, e1 = e1, e2 = e2)

    M_sb2_expected = [
      0.96940324     0.98982156     0.00000000     0.00000000     0.00000000     0.12337602  
      -0.06111586     0.96915934     0.00000000     0.00000000     0.00000000     0.24643492  
      0.00000000     0.00000000     1.00000000     1.00004154     0.00000000     0.00000000 
      0.00000000     0.00000000     0.00024545     1.00024546     0.00000000     0.00000000  
      -0.24643504    -0.12435558     0.00000000     0.00000000     1.00000000     0.48965743 
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000
    ]
    test_matrix(sb2, 1, M_sb2_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  =#

  # Linear Conbined Test
  L_b = 1.0
  e1 = 0.0
  e2 = 0.0
  K1 = 0.25

  #g = 0, gtot = 0, dg = 0, e = 0, k1 = 0.25,reduce to quad
    gtot = 0
    g = 0.0
    Kn1 = 0.25
    B0 = gtot * brho_ref
    Bn1 = Kn1 * brho_ref
    cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
    M_cb_expected = [
      0.87758256189037276     0.95885107720840601     0.00000000000000000     0.00000000000000000     0.00000000000000000     0.00000000000000000  
      -0.23971276930210150    0.87758256189037276     0.00000000000000000     0.00000000000000000     0.00000000000000000     0.00000000000000000
      0.00000000000000000     0.00000000000000000     1.12762596520638070     1.04219061098749480     0.00000000000000000     0.00000000000000000  
      0.00000000000000000     0.00000000000000000     0.26054765274687369     1.1276259652063807      0.00000000000000000     0.00000000000000000
      0.00000000000000000     0.00000000000000000     0.00000000000000000     0.00000000000000000     1.00000000     0.50000000  
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000 
    ]
    test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  #g = 0, e = 0, k1 = -0.25,should reduce to quad
    Kn1 = - 0.25
    Bn1 = Kn1 * brho_ref
    cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
    M_cb_expected = [
      1.1276259652063807     1.0421906109874948     0.00000000     0.00000000     0.00000000     0.00000000
      0.26054765274687369     1.1276259652063807     0.00000000     0.00000000     0.00000000     0.00000000
      0.00000000     0.00000000     0.87758256189037276     0.95885107720840601     0.00000000     0.00000000
      0.00000000     0.00000000    -0.23971276930210150     0.87758256189037276     0.00000000     0.00000000
      0.00000000     0.00000000     0.00000000     0.00000000     1.00000000     0.50000000
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000
    ]
    test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  # dg = 0, g = 0.25, k1 = 0.25, e = 0
    gtot = 0.25
    g = 0.25
    Kn1 = 0.25
    B0 = gtot * brho_ref
    Bn1 = Kn1 * brho_ref
    cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
    M_cb_expected = [
      8.4777686059853008e-01     9.4872443988117094e-01     0.0000000000000000     0.00000000     0.00000000     1.2177851152117594e-01 
      -2.9647638746286598e-01    8.4777686059853008e-01     0.0000000000000000     0.00000000     0.00000000     2.3718110997029276e-01
      0.00000000     0.00000000     1.1276259652063807      1.0421906109874948     0.00000000     0.00000000
      0.00000000     0.00000000     2.6054765274687369e-01  1.1276259652063807    0.00000000     0.00000000
      -2.3718110997029274e-01   -1.2177851152117591e-01     0.0000000000000000    0.00000000     1.00000000     4.8974488797623406e-01
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000
    ]
    test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  
  # dg = - 0.01, g = 0.25, k1 = 0.25, e = 0
    gtot = 0.24
    B0 = gtot * brho_ref
    cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
    M_cb_expected = [
      8.4896301853855183e-01     9.4912828113215142e-01     0.00000000     0.00000000     0.00000000     1.1705837590195874e-01  
      -2.9422976715096699e-01    8.4896301853855183e-01     0.00000000     0.00000000     0.00000000     2.3778289659600588e-01
      0.00000000     0.00000000  1.1276259652063807       1.0421906109874948     0.00000000     0.00000000
      0.00000000     0.00000000  2.6054765274687369e-01   1.1276259652063807    0.00000000     0.00000000
      -2.3631094433568914e-01    -1.2630823977784383e-01     0.00000000     0.00000000     1.00000000     4.8939468394148861e-01
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000
    ]
    test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)
    
  # dg = - 0.25, g = 0.25, k1 = 0.25, e = 0
    gtot = 0
    B0 = gtot * brho_ref
    cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
    M_cb_expected = [
      8.7758256189037276e-01     9.5885107720840601e-01     0.00000000     0.00000000     0.00000000     2.5610534585764899e-03 
      -2.3971276930210150e-01    8.7758256189037276e-01     0.00000000     0.00000000     0.00000000      2.4987133371685566e-01
      0.00000000     0.00000000    1.1276259652063807    1.0421906109874948     0.00000000     0.00000000 
      0.00000000     0.00000000    2.6054765274687369e-01     1.1276259652063807     0.00000000     0.00000000 
      -2.1989664240308857e-01    -2.3734186164259230e-01     0.00000000     0.00000000     1.00000000     4.9937479148421648e-01 
      0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000  
    ]
    test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)

    # dg = 0, g = 0.25, k1 = 0.25, e1 = e1 = 0.3
      gtot = 0.25
      g = 0.25
      B0 = gtot * brho_ref
      e1 = 0.3
      e2 = 0.3
      cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
      M_cb_expected = [
        9.2114557563498800e-01     9.4872443988117094e-01    0.00000000     0.00000000     0.00000000     1.2177851152117594e-01
        -1.5967842939416738e-01    9.2114557563498800e-01    0.00000000     0.00000000     0.00000000     2.4659873697954351e-01
        0.00000000     0.00000000     1.0470291314610722e+00     1.0421906109874948e+00     0.00000000     0.00000000
        0.00000000     0.00000000     9.2372739797483061e-02     1.0470291314610722e+00    0.00000000     0.00000000
        -2.4659873697954346e-01    -1.2177851152117591e-01    0.00000000     0.00000000     1.00000000     4.8974488797623406e-01
        0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000  
        ]
      test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)
    
    # # dg = -0.25 , g = 0.25, k1 = 0.25, e1 = e1 = 0.3
      #   g = 0.25
      #   gtot = 0.0
      #   Kn1 = 0.25
      #   B0 = gtot * brho_ref
      #   Bn1 = Kn1 * brho_ref
      #   cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
      #   M_cb_expected = [
          
      #   ]
      #   test_matrix(cb, 1, M_cb_expected, beta_gamma_ref=beta_gamma_ref, species=species)
      
    # dg = -0.15, g = 0.25  Kn1 =  0.25, e1 = 0.01, e2 = 0.9
        # g = 0.25
        # gtot = 0.10
        # Kn1 = 0.25
        # e1 = 0.01
        # e2 =0.9
        # B0 = gtot * brho_ref
        # Bn1 = Kn1*brho_ref
        # cb = Linear.Combined(L = L_b, B0 = B0, Bn1 = Bn1, g = g, e1 = e1, e2 = e2)
        # M_cb_expected = [
        #   0.8665771   0.9547928   0.0000000   0.0000000   0.0000000   0.0505521
        #   -0.1524999   0.9859413   0.0000000   0.0000000   0.0000000   0.2517563
        #    0.0000000   0.0000000   1.1265837   1.0421906   0.0000000   0.0000000
        #    0.0000000   0.0000000   0.1174526   0.9962935   0.0000000   0.0000000
        #   -0.2258755  -0.1905337   0.0000000   0.0000000   1.0000000   0.4908776
        #    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000
        # ] 

    end
#include("linear.jl")
include("element_tests.jl")