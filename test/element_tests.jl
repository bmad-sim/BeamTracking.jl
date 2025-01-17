using BeamTracking
#using AtomicAndPhysicalConstants
using Test

# define elements
# -- drifts
ld1 = 0.15;  # m
ld2 = 0.75;  # m
ld3 = 2.00;  # m
dr1 = MatrixKick.Drift(L = ld1)
dr2 = MatrixKick.Drift(L = ld2)
dr3 = MatrixKick.Drift(L = ld3)

# -- quadrupoles
lq1 = 0.05;  # m
lq2 = 0.30;  # m
lq3 = 0.75;  # m
lq4 = 1.20;  # m
gr1 = 0.20;  # T/m
gr2 = 0.60;  # T/m
gr3 = 1.20;  # T/m
gr4 = 3.50;  # T/m
qf1 = MatrixKick.Quadrupole(L = lq1, Bn1 = +gr1)
qd1 = MatrixKick.Quadrupole(L = lq1, Bn1 = -gr1)
qf2 = MatrixKick.Quadrupole(L = lq2, Bn1 = +gr2)
qd2 = MatrixKick.Quadrupole(L = lq2, Bn1 = -gr2)
qf3 = MatrixKick.Quadrupole(L = lq3, Bn1 = +gr3)
qd3 = MatrixKick.Quadrupole(L = lq3, Bn1 = -gr3)
qf4 = MatrixKick.Quadrupole(L = lq4, Bn1 = +gr4)
qd4 = MatrixKick.Quadrupole(L = lq4, Bn1 = -gr4)

# define beams
# -- species
e_minus = Species("electron")
p_plus =  Species("proton")
mec2 = massof(e_minus) # 0.51099895069 MeV
mpc2 = massof(p_plus)  # 938.27208943 MeV
# -- kinetic energy
ek1 =   5.e3;  # eV
ek2 =   1.e6;  # eV
ek3 =   1.e9;  # eV
ek4 = 250.e9;  # eV
# -- βγ
bg1 = sqrt(ek1 / mec2 * (ek1 / mec2 + 2))
bg2 = sqrt(ek2 / mec2 * (ek2 / mec2 + 2))
bg3 = sqrt(ek3 / mec2 * (ek3 / mec2 + 2))
bg4 = sqrt(ek4 / mpc2 * (ek4 / mpc2 + 2))
# -- pc
pc1 = mec2 * bg1
pc2 = mec2 * bg2
pc3 = mec2 * bg3
pc4 = mpc2 * bg4
# -- beams
# beam1.initial
xi  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ]
pxi = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ]
yi  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ]
pyi = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ]
zi  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ]
pzi = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ]

# beam1.dr1.final
xf_dr1  = [ 0.,  0.,                     0.,                     2.e-3,                  2.1123876489808660e-3, -2.1123876489808660e-3 ]
pxf_dr1 = [ 0.,  0.,                     0.,                     0.,                     7.5e-4,                -7.5e-4                ]
yf_dr1  = [ 0.,  0.,                     0.,                     1.e-3,                  1.0449550595923462e-3, -1.0449550595923462e-3 ]
pyf_dr1 = [ 0.,  0.,                     0.,                     0.,                     3.e-4,                 -3.e-4                 ]
zf_dr1  = [ 0.,  1.4710284465059844e-4, -1.4711135596840458e-4,  1.4710284465059844e-4,  1.4705400485512816e-4,  1.4705400485512816e-4 ]
pzf_dr1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.dr2.final
xf_dr2  = [ 0., 0., 0., 0.002, 0.0025619382449043287, -0.0025619382449043287 ]
pxf_dr2 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr2  = [ 0., 0., 0., 0.001, 0.0012247752979617315, -0.0012247752979617315 ]
pyf_dr2 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr2  = [ 0., 0.00008566359457101641, -0.00008589149602558208, 0.00008566359457101641, 0.00008541939559366522, 0.00008541939559366522 ]
pzf_dr2 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam3.dr3.final
xf_dr3  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
pxf_dr3 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr3  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
pyf_dr3 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr3  = [ 0., 5.209250185095532e-10, -5.224901403178192e-10, 5.209250185095532e-10, -6.506763479180273e-7, -6.506763479180273e-7 ]
pzf_dr3 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam4.dr4.final
xf_dr4  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
pxf_dr4 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr4  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
pyf_dr4 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr4  = [ 0., 2.7919184691863886e-8, -2.800306686850912e-8, 2.7919184691863886e-8, -6.232780882446728e-7, -6.232780882446728e-7 ]
pzf_dr4 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam1.qf1.final
xf_qf1  = [ 0.,  0.,                     0.,                     4.4834792779340600e-3,  4.5356143287504990e-3, -4.5356143287504990e-3 ]
pxf_qf1 = [ 0.,  0.,                     0.,                     1.1607821136240948e-1,  1.1776162121669208e-1, -1.1776162121669208e-1 ]
yf_qf1  = [ 0.,  0.,                     0.,                     1.2400905948673489e-4,  1.3427609296030678e-4, -1.3427609296030678e-4 ]
pyf_qf1 = [ 0.,  0.,                     0.,                    -2.8691666098954356e-2, -2.8653744321335432e-2,  2.8653744321335432e-2 ]
zf_qf1  = [ 0.,  4.903428155019947e-5,  -4.903711865613486e-5,  -4.8701323139842656e-5, -5.1605970700562340e-5, -5.1605970700562340e-5 ]
pzf_qf1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qd1.final
xf_qd1  = [ 0.,  0.,                     0.,                     2.4834727340278294e-4,  2.7409827330240064e-4, -2.7409827330240064e-4 ]
pxf_qd1 = [ 0.,  0.,                     0.,                    -5.7391734255615580e-2, -5.7299059109537503e-2,  5.7299059109537503e-2 ]
yf_qd1  = [ 0.,  0.,                     0.,                     2.2414071638535435e-3,  2.2621899982980444e-3, -2.2621899982980444e-3 ]
pyf_qd1 = [ 0.,  0.,                     0.,                     5.8033153149417160e-2,  5.8705242601074584e-2, -5.8705242601074584e-2 ]
zf_qd1  = [ 0.,  4.9034281550199470e-5, -4.9037118656134860e-5, -1.1242623339116514e-5, -1.1118475705649755e-5, -1.1118475705649755e-5 ]
pzf_qd1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.qf2.final
xf_qf2 =  [ 0.,  0.,                     0.,                     2.9271671401041872e-2,  3.0251479090608963e-2, -3.0251479090608963e-2 ]
pxf_qf2 = [ 0.,  0.,                     0.,                     3.2854870071095094e-1,  3.3959282727762820e-1, -3.3959282727762820e-1 ]
yf_qf2 =  [ 0.,  0.,                     0.,                    -9.7278287676729380e-4, -9.7883210731348450e-4,  9.7883210731348450e-4 ]
pyf_qf2 = [ 0.,  0.,                     0.,                     2.6415505265087600e-3,  2.3544440936341645e-3, -2.3544440936341645e-3 ]
zf_qf2 =  [ 0.,  3.4265437828406564e-5, -3.4356598410232830e-5, -2.3378443788359830e-3, -2.5012489357825860e-3, -2.5012489357825860e-3 ]
pzf_qf2 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.qd2.final
xf_qd2 =  [ 0.,  0.,                     0.,                    -0.0019464186216836795, -0.001961645985092862, 0.001961645985092862 ]
pxf_qd2 = [ 0.,  0.,                     0.,                     0.005200326165715375,   0.004472438532736971, -0.004472438532736971 ]
yf_qd2 =  [ 0.,  0.,                     0.,                     0.014608704236165976,   0.014997970393040412, -0.014997970393040412 ]
pyf_qd2 = [ 0.,  0.,                     0.,                     0.16398929979812307,    0.16837903611031713, -0.16837903611031713 ]
zf_qd2 =  [ 0.,  3.4265437828406564e-5, -3.435659841023283e-5, -5.899497183747271e-4, -6.222650472445767e-4, -6.222650472445767e-4 ]
pzf_qd2 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam3.qf3.final
xf_qf3 =  [ 0.,  0.,                      0.,                     0.002205479704035602,   0.0027865339903883585, -0.0027865339903883585 ]
pxf_qf3 = [ 0.,  0.,                      0.,                     0.0005576983180317967,  0.0013847532608658173, -0.0013847532608658173 ]
yf_qf3 =  [ 0.,  0.,                      0.,                     0.0009006624003320741,  0.0011179443238845692, -0.0011179443238845692 ]
pyf_qf3 = [ 0.,  0.,                      0.,                    -0.0002606852268330615,  9.513485195908057e-6, -9.513485195908057e-6 ]
zf_qf3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.630230762164262e-8,  -4.3665730716563075e-7, -4.3665730716563075e-7 ]
pzf_qf3 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam3.qd3.final
xf_qd3 =  [ 0.,  0.,                      0.,                     0.0018013248008431747,   0.0023445295165179687, -0.0023445295165179687 ]
pxf_qd3 = [ 0.,  0.,                      0.,                    -0.0005213704536906773,   0.0001541263391654702, -0.0001541263391654702 ]
yf_qd3 =  [ 0.,  0.,                      0.,                     0.0011027398519220506,   0.0013351614586315588, -0.0013351614586315588 ]
pyf_qd3 = [ 0.,  0.,                      0.,                     0.00027884915900320065,  0.0006096711218370137, -0.0006096711218370137 ]
zf_qd3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.4102076263621526e-8,  -1.6795543457606254e-7, -1.6795543457606254e-7 ]
pzf_qd3 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# beam4.qf4.final
xf_qf4 =  [ 0.,  0.,                      0.,                     0.0019939877700227713,   0.002892187842249561, -0.002892187842249561 ]
pxf_qf4 = [ 0.,  0.,                      0.,                    -0.000010025375230046909, 0.0007377200378071899, -0.0007377200378071899 ]
yf_qf4 =  [ 0.,  0.,                      0.,                     0.0010030091302526826,   0.0013630102695159176, -0.0013630102695159176 ]
pyf_qf4 = [ 0.,  0.,                      0.,                     5.022748546347732e-6,    0.0003059254879156335, -0.0003059254879156335 ]
zf_qf4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726401735327598e-8,  -3.698308917351584e-7, -3.698308917351584e-7 ]
pzf_qf4 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# beam4.qd4.final
xf_qd4 =  [ 0.,  0.,                      0.,                     0.002006018260505365,    0.00290602111426842, -0.00290602111426842 ]
pxf_qd4 = [ 0.,  0.,                      0.,                     0.000010045497092695465, 0.0007623023455299651, -0.0007623023455299651 ]
yf_qd4 =  [ 0.,  0.,                      0.,                     0.0009969938850113856,   0.0013562739161008426, -0.0013562739161008426 ]
pyf_qd4 = [ 0.,  0.,                      0.,                    -5.012687615023454e-6,    0.0002940854775943521, -0.0002940854775943521 ]
zf_qd4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726365460219405e-8,  -3.78176641902771e-7, -3.78176641902771e-7 ]
pzf_qd4 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# test individual elements
@testset "element_tests" begin
  # === drifts ===
  #
  # 5 keV electron
  beam1 = Bunch(species = e_minus, beta_gamma_ref = bg1,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam1, dr1);
  @test beam1.v.x  ≈  xf_dr1 (rtol=5.e-13)
  @test beam1.v.y  ≈  yf_dr1 (rtol=5.e-13)
  @test beam1.v.z  ≈  zf_dr1 (rtol=5.e-13)
  @test beam1.v.px == pxf_dr1
  @test beam1.v.py == pyf_dr1
  @test beam1.v.pz == pzf_dr1
  #
  # 1 MeV electron
  beam2 = Bunch(species = e_minus, beta_gamma_ref = bg2,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam2, dr2);
  @test beam2.v.x  ≈  xf_dr2 (rtol=5.e-13)
  @test beam2.v.y  ≈  yf_dr2 (rtol=5.e-13)
  @test beam2.v.z  ≈  zf_dr2 (rtol=5.e-13)
  @test beam2.v.px == pxf_dr2
  @test beam2.v.py == pyf_dr2
  @test beam2.v.pz == pzf_dr2
  #
  # 1 GeV electron
  beam3 = Bunch(species = e_minus, beta_gamma_ref = bg3,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam3, dr3);
  @test beam3.v.x  ≈  xf_dr3 (rtol=5.e-13)
  @test beam3.v.y  ≈  yf_dr3 (rtol=5.e-13)
  @test beam3.v.z  ≈  zf_dr3 (rtol=5.e-13)
  @test beam3.v.px == pxf_dr3
  @test beam3.v.py == pyf_dr3
  @test beam3.v.pz == pzf_dr3
  #
  # 250 GeV proton
  beam4 = Bunch(species = p_plus, beta_gamma_ref = bg4,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam4, dr3);
  @test beam4.v.x  ≈  xf_dr4 (rtol=5.e-13)
  @test beam4.v.y  ≈  yf_dr4 (rtol=5.e-13)
  @test beam4.v.z  ≈  zf_dr4 (rtol=5.e-13)
  @test beam4.v.px == pxf_dr4
  @test beam4.v.py == pyf_dr4
  @test beam4.v.pz == pzf_dr4

  # === quadrupoles ===
  #
  # 5 keV electron
  beam1 = Bunch(species = e_minus, beta_gamma_ref = bg1,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam1, qf1);
  @test beam1.v.x  ≈  xf_qf1  (rtol=5.e-13)
  @test beam1.v.px ≈  pxf_qf1 (rtol=5.e-13)
  @test beam1.v.y  ≈  yf_qf1  (rtol=5.e-13)
  @test beam1.v.py ≈  pyf_qf1 (rtol=5.e-13)
  @test beam1.v.z  ≈  zf_qf1  (rtol=5.e-13)
  @test beam1.v.pz == pzf_qf1
  beam1 = Bunch(species = e_minus, beta_gamma_ref = bg1,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam1, qd1);
  @test beam1.v.x  ≈  xf_qd1  (rtol=5.e-13)
  @test beam1.v.px ≈  pxf_qd1 (rtol=5.e-13)
  @test beam1.v.y  ≈  yf_qd1  (rtol=5.e-13)
  @test beam1.v.py ≈  pyf_qd1 (rtol=5.e-13)
  @test beam1.v.z  ≈  zf_qd1  (rtol=5.e-13)
  @test beam1.v.pz == pzf_qd1
  #
  # 1 MeV electron
  beam2 = Bunch(species = e_minus, beta_gamma_ref = bg2,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam2, qf2);
  @test beam2.v.x  ≈  xf_qf2  (rtol=5.e-13)
  @test beam2.v.px ≈  pxf_qf2 (rtol=5.e-13)
  @test beam2.v.y  ≈  yf_qf2  (rtol=5.e-13)
  @test beam2.v.py ≈  pyf_qf2 (rtol=5.e-13)
  @test beam2.v.z  ≈  zf_qf2  (rtol=5.e-13)
  @test beam2.v.pz == pzf_qf2
  beam2 = Bunch(species = e_minus, beta_gamma_ref = bg2,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam2, qd2);
  @test beam2.v.x  ≈  xf_qd2  (rtol=5.e-13)
  @test beam2.v.px ≈  pxf_qd2 (rtol=5.e-13)
  @test beam2.v.y  ≈  yf_qd2  (rtol=5.e-13)
  @test beam2.v.py ≈  pyf_qd2 (rtol=5.e-13)
  @test beam2.v.z  ≈  zf_qd2  (rtol=5.e-13)
  @test beam2.v.pz == pzf_qd2
  #
  # 1 GeV electron
  beam3 = Bunch(species = e_minus, beta_gamma_ref = bg3,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam3, qf3);
  @test beam3.v.x  ≈  xf_qf3  (rtol=5.e-13)
  @test beam3.v.px ≈  pxf_qf3 (rtol=5.e-13)
  @test beam3.v.y  ≈  yf_qf3  (rtol=5.e-13)
  @test beam3.v.py ≈  pyf_qf3 (rtol=5.e-13)
  @test beam3.v.z  ≈  zf_qf3  (rtol=5.e-13)
  @test beam3.v.pz == pzf_qf3
  beam3 = Bunch(species = e_minus, beta_gamma_ref = bg3,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam3, qd3);
  @test beam3.v.x  ≈  xf_qd3  (rtol=5.e-13)
  @test beam3.v.px ≈  pxf_qd3 (rtol=5.e-13)
  @test beam3.v.y  ≈  yf_qd3  (rtol=5.e-13)
  @test beam3.v.py ≈  pyf_qd3 (rtol=5.e-13)
  @test beam3.v.z  ≈  zf_qd3  (rtol=5.e-13)
  @test beam3.v.pz == pzf_qd3
  #
  # 250 GeV proton
  beam4 = Bunch(species = p_plus, beta_gamma_ref = bg4,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam4, qf4);
  @test beam4.v.x  ≈  xf_qf4  (rtol=5.e-13)
  @test beam4.v.px ≈  pxf_qf4 (rtol=5.e-13)
  @test beam4.v.y  ≈  yf_qf4  (rtol=5.e-13)
  @test beam4.v.py ≈  pyf_qf4 (rtol=5.e-13)
  @test beam4.v.z  ≈  zf_qf4  (rtol=5.e-13)
  @test beam4.v.pz == pzf_qf4
  beam4 = Bunch(species = p_plus, beta_gamma_ref = bg4,
               x = copy(xi), px = copy(pxi), y = copy(yi), py = copy(pyi), z = copy(zi), pz = copy(pzi))
  track!(beam4, qd4);
  @test beam4.v.x  ≈  xf_qd4  (rtol=5.e-13)
  @test beam4.v.px ≈  pxf_qd4 (rtol=5.e-13)
  @test beam4.v.y  ≈  yf_qd4  (rtol=5.e-13)
  @test beam4.v.py ≈  pyf_qd4 (rtol=5.e-13)
  @test beam4.v.z  ≈  zf_qd4  (rtol=5.e-13)
  @test beam4.v.pz == pzf_qd4


end

