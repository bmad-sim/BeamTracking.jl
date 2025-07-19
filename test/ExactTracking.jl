# test individual elements
@testset "ExactTracking" begin
    
  @testset "Particles" begin
    include("extras/tracking_values.jl")
    # ===  D R I F T  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_drift!, (β1, γsq1, 1/βγ1, ld1)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_drift!, (β2, γsq2, 1/βγ2, ld2)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_drift!, (β3, γsq3, 1/βγ3, ld3)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_drift!, (β4, γsq4, 1/βγ4, ld4)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr4 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr4 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr4 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #=
    # ===  Q U A D R U P O L E  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1,  gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1, -gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2,  gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2, -gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3,  gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3, -gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4,  gr4 / Bρ4, lq4)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4, -gr4 / Bρ4, lq4)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    =#
    # ===  T H I N - L E N S   K I C K  ===
    #
    # 5 keV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_k1 * cos(ra1) / Bρ1
    ks1 = bv_k1 * sin(ra1) / Bρ1
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k1,  kn1 * lk1,  ks1 * lk1, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k1, -kn1 * lk1, -ks1 * lk1, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_dk1 * cos(ra1) / Bρ1
    ks1 = bv_dk1 * sin(ra1) / Bρ1
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk1,  kn1 * lk1,  ks1 * lk1, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk1, -kn1 * lk1, -ks1 * lk1, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 1 MeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn2 = bv_k2 * cos(ra2) / Bρ2
    ks2 = bv_k2 * sin(ra2) / Bρ2
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k2,  kn2 * lk2,  ks2 * lk2, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k2, -kn2 * lk2, -ks2 * lk2, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn2 = bv_dk2 * cos(ra2) / Bρ2
    ks2 = bv_dk2 * sin(ra2) / Bρ2
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk2,  kn2 * lk2,  ks2 * lk2, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk2, -kn2 * lk2, -ks2 * lk2, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 1 GeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn3 = bv_k3 * cos(ra3) / Bρ3
    ks3 = bv_k3 * sin(ra3) / Bρ3
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k3,  kn3 * lk3,  ks3 * lk3, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k3, -kn3 * lk3, -ks3 * lk3, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn3 = bv_dk3 * cos(ra3) / Bρ3
    ks3 = bv_dk3 * sin(ra3) / Bρ3
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk3,  kn3 * lk3,  ks3 * lk3, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk3, -kn3 * lk3, -ks3 * lk3, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 250 GeV proton
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn4 = bv_k4 * cos(ra4) / Bρ4
    ks4 = bv_k4 * sin(ra4) / Bρ4
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k4,  kn4 * lk4,  ks4 * lk4, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k4, -kn4 * lk4, -ks4 * lk4, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn4 = bv_dk4 * cos(ra4) / Bρ4
    ks4 = bv_dk4 * sin(ra4) / Bρ4
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk4,  kn4 * lk4,  ks4 * lk4, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk4, -kn4 * lk4, -ks4 * lk4, 1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2

    # ===  M U L T I P O L E  ===
    #
    # 5 keV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_m1 * cos(ra1) / Bρ1
    ks1 = bv_m1 * sin(ra1) / Bρ1
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β1, γsq1, 1/βγ1, ms_m1,  kn1,  ks1, lm1)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β1, γsq1, 1/βγ1, ms_m1, -kn1, -ks1, lm1)))
    @test v[:,BeamTracking.XI]  ≈  xf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 1 MeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn2 = bv_m2 * cos(ra2) / Bρ2
    ks2 = bv_m2 * sin(ra2) / Bρ2
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β2, γsq2, 1/βγ2, ms_m2,  kn2,  ks2, lm2)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β2, γsq2, 1/βγ2, ms_m2, -kn2, -ks2, lm2)))
    @test v[:,BeamTracking.XI]  ≈  xf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 1 GeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn3 = bv_m3 * cos(ra3) / Bρ3
    ks3 = bv_m3 * sin(ra3) / Bρ3
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β3, γsq3, 1/βγ3, ms_m3,  kn3,  ks3, lm3)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β3, γsq3, 1/βγ3, ms_m3, -kn3, -ks3, lm3)))
    @test v[:,BeamTracking.XI]  ≈  xf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    #
    # 250 GeV proton
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn4 = bv_m4 * cos(ra4) / Bρ4
    ks4 = bv_m4 * sin(ra4) / Bρ4
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β4, γsq4, 1/βγ4, ms_m4,  kn4,  ks4, lm4)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(IntegrationTracking.dkd_multipole!, (β4, γsq4, 1/βγ4, ms_m4, -kn4, -ks4, lm4)))
    @test v[:,BeamTracking.XI]  ≈  xf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2

    #=
    # ===  S B E N D  ===
    #
    # 5 keV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_sbend!, (β1, Bρ1, hc1, b_1, ee1, ex1, la1)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb1
    #
    # 1 MeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_sbend!, (β2, Bρ2, hc2, b_2, ee2, ex2, la2)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb2
    #
    # 1 GeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_sbend!, (β3, Bρ3, hc3, b_3, ee3, ex3, la3)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb3
    #
    # 250 GeV proton
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.exact_sbend!, (β4, Bρ4, hc4, b_4, ee4, ex4, la4)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb4
    =#
  end
    
  @testset "Utility functions" begin
    dx_rot = -0.1
    dy_rot = -0.1
    dz_rot = 0.2

    W = [cos(dy_rot) 0 sin(dy_rot); 0 1 0; -sin(dy_rot) 0 cos(dy_rot)] *
        [1 0 0; 0 cos(dx_rot) -sin(dx_rot); 0 sin(dx_rot) cos(dx_rot)] *
        [cos(dz_rot) -sin(dz_rot) 0; sin(dz_rot) cos(dz_rot) 0; 0 0 1]

    # Test w_matrix function
    @test all(ExactTracking.w_matrix(dx_rot, dy_rot, dz_rot) .== W)

    Winv = [cos(dz_rot) sin(dz_rot) 0; -sin(dz_rot) cos(dz_rot) 0; 0 0 1] *
            [1 0 0; 0 cos(dx_rot) sin(dx_rot); 0 -sin(dx_rot) cos(dx_rot)] *
            [cos(dy_rot) 0 -sin(dy_rot); 0 1 0; sin(dy_rot) 0 cos(dy_rot)]

    # Test w_inv_matrix function
    @test all(ExactTracking.w_inv_matrix(dx_rot, dy_rot, dz_rot) .== Winv)
  end

  @testset "Kernels" begin
    function patch_args(::Type{T}) where {T}
        p0c = T(10e6)
        mc2 = T(ELECTRON.mass)
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        dt = T(1e-9)
        dx = T(2)
        dy = T(3)
        dz = T(4)
        winv = ExactTracking.w_inv_matrix(T(-5),T(6),T(7))
        L = winv[3,1]*dx + winv[3,2]*dy + winv[3,3]*dz
        return beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, winv, L
    end

    function patch_norot_args(::Type{T}) where {T}
        p0c = T(10e6)
        mc2 = T(ELECTRON.mass)
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        dt = T(4e-9)
        dx = T(1)
        dy = T(2)
        dz = T(3)
        L = dz
        return beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, nothing, L
    end

    function drift_args(::Type{T}) where {T}
        L = T(1)
        p0c = T(10e6)
        mc2 = T(ELECTRON.mass)
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        return beta_0, gamsqr_0, tilde_m, L
    end
    
    function solenoid_args(::Type{T}) where {T}
        L = T(1)
        ks = T(2)
        p0c = T(10e6)
        mc2 = T(ELECTRON.mass)
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        return ks, beta_0, gamsqr_0, tilde_m, L
    end

    # Scalar parameters
    test_map("bmad_maps/patch.jl",       KernelCall(ExactTracking.patch!, patch_args(Float64));                           tol=5e-10)
    test_map("bmad_maps/patch_norot.jl", KernelCall(ExactTracking.patch!, patch_norot_args(Float64));                     tol=1e-9 )
    test_map("bmad_maps/drift.jl",       KernelCall(ExactTracking.ExactTracking.exact_drift!, drift_args(Float64));       tol=5e-10)
    test_map("bmad_maps/solenoid.jl",    KernelCall(ExactTracking.ExactTracking.exact_solenoid!, solenoid_args(Float64)); tol=5e-10)

    # GTPSA parameters
    test_map("bmad_maps/patch.jl",       KernelCall(ExactTracking.patch!, patch_args(TPS64{D10}));                           tol=5e-10)
    test_map("bmad_maps/patch_norot.jl", KernelCall(ExactTracking.patch!, patch_norot_args(TPS64{D10}));                     tol=1e-9 )
    test_map("bmad_maps/drift.jl",       KernelCall(ExactTracking.ExactTracking.exact_drift!, drift_args(TPS64{D10}));       tol=5e-10)
    test_map("bmad_maps/solenoid.jl",    KernelCall(ExactTracking.ExactTracking.exact_solenoid!, solenoid_args(TPS64{D10})); tol=5e-10)

    # Update P0 for a bunch
    v0 = [0.0 0.1 0.0 0.1 0.0 0.0]
    b0 = BunchView([State.Alive], copy(v0), nothing)
    kernel_call = KernelCall(ExactTracking.update_P0!, (1.0, 2.0))
    BeamTracking.launch!(b0, kernel_call)
    @test b0.v == [0.0 0.05 0.0 0.05 0.0 -0.5]
    @test_opt kernel_call.kernel(1, b0, kernel_call.args...)
    @test @ballocated(BeamTracking.launch!(b0, $kernel_call; use_KA=false), 
    setup=(b0 = BunchView([State.Alive], copy($v0), nothing))) == 0
  end
end