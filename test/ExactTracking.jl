c_light = 299792458.0
degree = π / 180


# ========== define element paameters ==========

# -- drifts
ld1 = 0.15  # m
ld2 = 0.75  # m
ld3 = 2.00  # m
ld4 = 2.00  # m

# -- quadrupoles
lq1 = 0.05  # m
lq2 = 0.30  # m
lq3 = 0.75  # m
lq4 = 1.20  # m

gr1 = 0.20  # T/m
gr2 = 0.60  # T/m
gr3 = 1.20  # T/m
gr4 = 3.50  # T/m

# -- thin-lens kicks
lk1 = 0.05  # m
lk2 = 0.30  # m
lk3 = 0.75  # m
lk4 = 2.25  # m

fact0_5 = [ 1, 1, 2, 6, 24, 120 ]
fact1_5 = [    1, 2, 6, 24, 120 ]
fact2_5 = [       2, 6, 24, 120 ]

ra1 = 15degree
ms_k1  = [          2,    3,    4,    5,    6   ]
bv_k1  = [         0.20, 0.18, 0.16, 0.14, 0.12 ] .* fact1_5  # T/m^{n-1}
ms_dk1 = [  1,      2,    3,    4,    5,    6   ]
bv_dk1 = [ 0.4e-3, 0.20, 0.18, 0.16, 0.14, 0.12 ] .* fact0_5  # T/m^{n-1}

ra2 = 12degree
ms_k2  = [          2,    3,    4,    5,    6   ]
bv_k2  = [         0.60, 0.40, 0.30, 0.20, 0.10 ] .* fact1_5  # T·m^{n-1}
ms_dk2 = [  1,      2,    3,    4,    5,    6   ]
bv_dk2 = [ 1.2e-3, 0.60, 0.40, 0.30, 0.20, 0.10 ] .* fact0_5  # T/m^{n-1}

ra3 = 9degree
ms_k3  = [          2,    3,    4,    5,    6   ]
bv_k3  = [         1.20, 1.00, 0.80, 0.60, 0.40 ] .* fact1_5  # T·m^{n-1}
ms_dk3 = [  1,      2,    3,    4,    5,    6   ]
bv_dk3 = [ 0.03,   1.20, 1.00, 0.80, 0.60, 0.40 ] .* fact0_5  # T/m^{n-1}

ra4 = 6degree
ms_k4  = [          2,    3,    4,    5,    6   ]
bv_k4  = [         3.60, 3.00, 2.40, 1.80, 1.20 ] .* fact1_5  # T·m^{n-1}
ms_dk4 = [  1,      2,    3,    4,    5,    6   ]
bv_dk4 = [ 0.39,   3.60, 3.00, 2.40, 1.80, 1.20 ] .* fact0_5  # T/m^{n-1}

# -- multipoles
lm1 = 0.05  # m
lm2 = 0.30  # m
lm3 = 0.75  # m
lm4 = 2.25  # m

ra1 = 15degree
ms_m1  = [                3,    4,    5,    6   ]
bv_m1  = [               0.18, 0.16, 0.14, 0.12 ] .* fact2_5  # T/m^{n-1}

ra2 = 12degree
ms_m2  = [                3,    4,    5,    6   ]
bv_m2  = [               0.40, 0.30, 0.20, 0.10 ] .* fact2_5  # T/m^{n-1}

ra3 = 9degree
ms_m3  = [                3,    4,    5,    6   ]
bv_m3  = [               1.00, 0.80, 0.60, 0.40 ] .* fact2_5  # T/m^{n-1}

ra4 = 6degree
ms_m4  = [                3,    4,    5,    6   ]
bv_m4  = [               3.00, 2.40, 1.80, 1.20 ] .* fact2_5  # T/m^{n-1}


# ========== define beams ==========

# -- species
#e_minus = Species("electron")
#p_plus =  Species("proton")
#mec2 = BeamTracking.massof(e_minus) # 0.51099895069 MeV
#mpc2 = BeamTracking.massof(p_plus)  # 938.27208943 MeV
mec2 = 0.51099895069e6 # eV
mpc2 = 938.27208943e6  # eV

# -- kinetic energy
ek1 =   5.e3;  # eV
ek2 =   1.e6;  # eV
ek3 =   1.e9;  # eV
ek4 = 250.e9;  # eV

# -- βγ
βγ1 = sqrt(ek1 / mec2 * (ek1 / mec2 + 2))
βγ2 = sqrt(ek2 / mec2 * (ek2 / mec2 + 2))
βγ3 = sqrt(ek3 / mec2 * (ek3 / mec2 + 2))
βγ4 = sqrt(ek4 / mpc2 * (ek4 / mpc2 + 2))

# -- γ^2
γsq1 = 1.0 + βγ1^2
γsq2 = 1.0 + βγ2^2
γsq3 = 1.0 + βγ3^2
γsq4 = 1.0 + βγ4^2

# -- β
β1 = βγ1 / sqrt(γsq1)
β2 = βγ2 / sqrt(γsq2)
β3 = βγ3 / sqrt(γsq3)
β4 = βγ4 / sqrt(γsq4)

# -- pc
pc1 = mec2 * βγ1
pc2 = mec2 * βγ2
pc3 = mec2 * βγ3
pc4 = mpc2 * βγ4

# -- Bρ
Bρ1 = -pc1 / c_light
Bρ2 = -pc2 / c_light
Bρ3 = -pc3 / c_light
Bρ4 = +pc4 / c_light


# ========== specify initial conditions ==========

# beam.initial
xi  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ]
pxi = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ]
yi  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ]
pyi = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ]
zi  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ]
pzi = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ]

# beam2.initial
xi2  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200,  0.00200, -0.00200 ]
pxi2 = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075,  0.00075, -0.00075 ]
yi2  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100,  0.00100, -0.00100 ]
pyi2 = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030,  0.00030, -0.00030 ]
zi2  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000,  0.00000,  0.00000 ]
pzi2 = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100, -0.00100, -0.00100 ]


# ========== specify final conditions ==========

# -- drifts

# beam.dr1.final
xf_dr1  = [ 0.,  0.,                     0.,                     2.e-3,                  2.1123876489808660e-3, -2.1123876489808660e-3 ]
yf_dr1  = [ 0.,  0.,                     0.,                     1.e-3,                  1.0449550595923462e-3, -1.0449550595923462e-3 ]
zf_dr1  = [ 0.,  1.4710284465059844e-4, -1.4711135596840458e-4,  1.4710284465059844e-4,  1.4705400485512816e-4,  1.4705400485512816e-4 ]

# beam.dr2.final
xf_dr2  = [ 0., 0., 0., 0.002, 0.0025619382449043287, -0.0025619382449043287 ]
yf_dr2  = [ 0., 0., 0., 0.001, 0.0012247752979617315, -0.0012247752979617315 ]
zf_dr2  = [ 0., 0.00008566359457101641, -0.00008589149602558208, 0.00008566359457101641, 0.00008541939559366522, 0.00008541939559366522 ]

# beam.dr3.final
xf_dr3  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
yf_dr3  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
zf_dr3  = [ 0., 5.209250185095532e-10, -5.224901403178192e-10, 5.209250185095532e-10, -6.506763479180273e-7, -6.506763479180273e-7 ]

# beam.dr4.final
xf_dr4  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
yf_dr4  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
zf_dr4  = [ 0., 2.7919184691863886e-8, -2.800306686850912e-8, 2.7919184691863886e-8, -6.232780882446728e-7, -6.232780882446728e-7 ]

# -- quadrupoles

# beam.qf1.final
xf_qf1  = [ 0.,  0.,                      0.,                     4.4834792779340600e-3,  4.5356143287504990e-3, -4.5356143287504990e-3 ]
pxf_qf1 = [ 0.,  0.,                      0.,                     1.1607821136240948e-1,  1.1776162121669208e-1, -1.1776162121669208e-1 ]
yf_qf1  = [ 0.,  0.,                      0.,                     1.2400905948673489e-4,  1.3427609296030678e-4, -1.3427609296030678e-4 ]
pyf_qf1 = [ 0.,  0.,                      0.,                    -2.8691666098954356e-2, -2.8653744321335432e-2,  2.8653744321335432e-2 ]
zf_qf1  = [ 0.,  4.903428155019947e-5,   -4.903711865613486e-5,  -4.8701323139842656e-5, -5.1605970700562340e-5, -5.1605970700562340e-5 ]

# beam.qd1.final
xf_qd1  = [ 0.,  0.,                      0.,                     2.4834727340278294e-4,  2.7409827330240064e-4, -2.7409827330240064e-4 ]
pxf_qd1 = [ 0.,  0.,                      0.,                    -5.7391734255615580e-2, -5.7299059109537503e-2,  5.7299059109537503e-2 ]
yf_qd1  = [ 0.,  0.,                      0.,                     2.2414071638535435e-3,  2.2621899982980444e-3, -2.2621899982980444e-3 ]
pyf_qd1 = [ 0.,  0.,                      0.,                     5.8033153149417160e-2,  5.8705242601074584e-2, -5.8705242601074584e-2 ]
zf_qd1  = [ 0.,  4.9034281550199470e-5,  -4.9037118656134860e-5, -1.1242623339116514e-5, -1.1118475705649755e-5, -1.1118475705649755e-5 ]

# beam.qf2.final
xf_qf2 =  [ 0.,  0.,                      0.,                     2.9271671401041872e-2,  3.0251479090608963e-2, -3.0251479090608963e-2 ]
pxf_qf2 = [ 0.,  0.,                      0.,                     3.2854870071095094e-1,  3.3959282727762820e-1, -3.3959282727762820e-1 ]
yf_qf2 =  [ 0.,  0.,                      0.,                    -9.7278287676729380e-4, -9.7883210731348450e-4,  9.7883210731348450e-4 ]
pyf_qf2 = [ 0.,  0.,                      0.,                     2.6415505265087600e-3,  2.3544440936341645e-3, -2.3544440936341645e-3 ]
zf_qf2 =  [ 0.,  3.4265437828406564e-5,  -3.4356598410232830e-5, -2.3378443788359830e-3, -2.5012489357825860e-3, -2.5012489357825860e-3 ]

# beam.qd2.final
xf_qd2 =  [ 0.,  0.,                      0.,                    -0.0019464186216836795, -0.001961645985092862, 0.001961645985092862 ]
pxf_qd2 = [ 0.,  0.,                      0.,                     0.005200326165715375,   0.004472438532736971, -0.004472438532736971 ]
yf_qd2 =  [ 0.,  0.,                      0.,                     0.014608704236165976,   0.014997970393040412, -0.014997970393040412 ]
pyf_qd2 = [ 0.,  0.,                      0.,                     0.16398929979812307,    0.16837903611031713, -0.16837903611031713 ]
zf_qd2 =  [ 0.,  3.4265437828406564e-5,  -3.435659841023283e-5,  -5.899497183747271e-4,  -6.222650472445767e-4, -6.222650472445767e-4 ]

# beam.qf3.final
xf_qf3 =  [ 0.,  0.,                      0.,                     0.002205479704035602,   0.0027865339903883585, -0.0027865339903883585 ]
pxf_qf3 = [ 0.,  0.,                      0.,                     0.0005576983180317967,  0.0013847532608658173, -0.0013847532608658173 ]
yf_qf3 =  [ 0.,  0.,                      0.,                     0.0009006624003320741,  0.0011179443238845692, -0.0011179443238845692 ]
pyf_qf3 = [ 0.,  0.,                      0.,                    -0.0002606852268330615,  9.513485195908057e-6, -9.513485195908057e-6 ]
zf_qf3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.630230762164262e-8,  -4.3665730716563075e-7, -4.3665730716563075e-7 ]

# beam.qd3.final
xf_qd3 =  [ 0.,  0.,                      0.,                     0.0018013248008431747,   0.0023445295165179687, -0.0023445295165179687 ]
pxf_qd3 = [ 0.,  0.,                      0.,                    -0.0005213704536906773,   0.0001541263391654702, -0.0001541263391654702 ]
yf_qd3 =  [ 0.,  0.,                      0.,                     0.0011027398519220506,   0.0013351614586315588, -0.0013351614586315588 ]
pyf_qd3 = [ 0.,  0.,                      0.,                     0.00027884915900320065,  0.0006096711218370137, -0.0006096711218370137 ]
zf_qd3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.4102076263621526e-8,  -1.6795543457606254e-7, -1.6795543457606254e-7 ]

# beam.qf4.final
xf_qf4 =  [ 0.,  0.,                      0.,                     0.0019939877700227713,   0.002892187842249561, -0.002892187842249561 ]
pxf_qf4 = [ 0.,  0.,                      0.,                    -0.000010025375230046909, 0.0007377200378071899, -0.0007377200378071899 ]
yf_qf4 =  [ 0.,  0.,                      0.,                     0.0010030091302526826,   0.0013630102695159176, -0.0013630102695159176 ]
pyf_qf4 = [ 0.,  0.,                      0.,                     5.022748546347732e-6,    0.0003059254879156335, -0.0003059254879156335 ]
zf_qf4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726401735327598e-8,  -3.698308917351584e-7, -3.698308917351584e-7 ]

# beam.qd4.final
xf_qd4 =  [ 0.,  0.,                      0.,                     0.002006018260505365,    0.00290602111426842, -0.00290602111426842 ]
pxf_qd4 = [ 0.,  0.,                      0.,                     0.000010045497092695465, 0.0007623023455299651, -0.0007623023455299651 ]
yf_qd4 =  [ 0.,  0.,                      0.,                     0.0009969938850113856,   0.0013562739161008426, -0.0013562739161008426 ]
pyf_qd4 = [ 0.,  0.,                      0.,                    -5.012687615023454e-6,    0.0002940854775943521, -0.0002940854775943521 ]
zf_qd4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726365460219405e-8,  -3.78176641902771e-7, -3.78176641902771e-7 ]

# -- thin-lens kicks

# beam2.kp1.final
pxf_kp1 = [  0.,  0.,     0.,     0.07006321275745035,  0.07081321275745035, -0.07067295798136997,  0.07081321275745035, -0.07067295798136997 ]
pyf_kp1 = [  0.,  0.,     0.,    -0.06224157687021319, -0.06194157687021319,  0.06159214864638122, -0.06194157687021319,  0.06159214864638122 ]

# beam2.kn1.final
pxf_kn1 = [  0.,  0.,     0.,    -0.07006321275745035, -0.06931321275745035,  0.06917295798136996, -0.06931321275745035,  0.06917295798136996 ]
pyf_kn1 = [  0.,  0.,     0.,     0.06224157687021319,  0.0625415768702132,  -0.06219214864638122,  0.0625415768702132,  -0.06219214864638122 ]

# beam2.dkp1.final
pxf_dkp1 = [  0.08082108872791081,   0.08082108872791081,   0.08082108872791081,   0.15088430148536117,  0.15163430148536117,  0.01014813074654084,   0.15163430148536117,  0.01014813074654084  ]
pyf_dkp1 = [ -0.021655945456047817, -0.021655945456047817, -0.021655945456047817, -0.083897522326261,   -0.08359752232626101,  0.039936203190333405, -0.08359752232626101,  0.039936203190333405 ]

# beam2.dkn1.final
pxf_dkn1 = [ -0.08082108872791081,  -0.08082108872791081,  -0.08082108872791081,  -0.15088430148536117, -0.15013430148536117, -0.011648130746540841, -0.15013430148536117, -0.011648130746540841 ]
pyf_dkn1 = [  0.021655945456047817,  0.021655945456047817,  0.021655945456047817,  0.083897522326261,    0.084197522326261,   -0.04053620319033341,   0.084197522326261,   -0.04053620319033341  ]

# beam2.kp2.final
pxf_kp2 = [  0.,  0.,     0.,     0.06640298429951919,  0.06715298429951919, -0.06704658530538163,  0.06715298429951919, -0.06704658530538163 ]
pyf_kp2 = [  0.,  0.,     0.,    -0.05301509095771739, -0.05271509095771739,  0.05248555704730759, -0.05271509095771739,  0.05248555704730759 ]

# beam2.kn2.final
pxf_kn2 = [  0.,  0.,     0.,    -0.06640298429951919, -0.06565298429951918,   0.065546585305381630, -0.06565298429951918,   0.06554658530538163  ]
pyf_kn2 = [  0.,  0.,     0.,     0.05301509095771739,  0.053315090957717394, -0.053085557047307594,  0.053315090957717394, -0.053085557047307594 ]

# beam2.dkp2.final
pxf_dkp2 = [  0.07423987764364487,  0.07423987764364487,  0.07423987764364487,  0.14064286194316405,  0.14139286194316406,  0.0071932923382632364,  0.14139286194316406,  0.0071932923382632364 ]
pyf_dkp2 = [ -0.01578017313073629, -0.01578017313073629, -0.01578017313073629, -0.06879526408845368, -0.06849526408845369,  0.0367053839165713,    -0.06849526408845369,  0.0367053839165713    ]

# beam2.dkn2.final
pxf_dkn2 = [ -0.07423987764364487, -0.07423987764364487, -0.07423987764364487, -0.14064286194316405, -0.13989286194316405, -0.008693292338263237,  -0.13989286194316405, -0.008693292338263237 ]
pyf_dkn2 = [  0.01578017313073629,  0.01578017313073629,  0.01578017313073629,  0.06879526408845368,  0.06909526408845368, -0.0373053839165713,     0.06909526408845368, -0.0373053839165713   ]

# beam2.kp3.final
pxf_kp3 = [  0.,  0.,     0.,     0.0004910493596606422,  0.0012410493596606421, -0.0012399988296393268,    0.0012410493596606421, -0.0012399988296393268 ]
pyf_kp3 = [  0.,  0.,     0.,    -0.0003517236886733822, -0.0000517236886733822,  0.000049737047965373754, -0.0000517236886733822,  0.000049737047965373754 ]

# beam2.kn3.final
pxf_kn3 = [  0.,  0.,     0.,    -0.0004910493596606422,  0.0002589506403393578, -0.0002600011703606731,  0.0002589506403393578, -0.0002600011703606731 ]
pyf_kn3 = [  0.,  0.,     0.,     0.0003517236886733822,  0.0006517236886733821, -0.0006497370479653737,  0.0006517236886733821, -0.0006497370479653737 ]

# beam2.dkp3.final
pxf_dkp3 = [  0.006658882282791728,   0.006658882282791728,   0.006658882282791728,   0.007149931642452371,   0.007899931642452371,   0.005418883453152402,   0.007899931642452371,   0.005418883453152402  ]
pyf_dkp3 = [ -0.0010546633435469385, -0.0010546633435469385, -0.0010546633435469385, -0.0014063870322203209, -0.001106387032220321,  -0.0010049262955815648, -0.001106387032220321,  -0.0010049262955815648 ]

# beam2.dkn3.final
pxf_dkn3 = [ -0.006658882282791728,  -0.006658882282791728,  -0.006658882282791728,  -0.007149931642452371,  -0.006399931642452371,  -0.006918883453152402,  -0.006399931642452371,  -0.006918883453152402  ]
pyf_dkn3 = [  0.0010546633435469385,  0.0010546633435469385,  0.0010546633435469385,  0.0014063870322203209,  0.0017063870322203208,  0.0004049262955815648,  0.0017063870322203208,  0.0004049262955815648 ]

# beam2.kp4.final
pxf_kp4 = [  0.,  0.,  0., -0.000018257195229563586,  0.0007317428047704365,  -0.0007317841812471834,   0.0007317428047704365,  -0.0007317841812471834  ]
pyf_kp4 = [  0.,  0.,  0.,  0.00001168174301104629,   0.00031168174301104626, -0.00031161252507745754,  0.00031168174301104626, -0.00031161252507745754 ]

# beam2.kn4.final
pxf_kn4 = [  0.,  0.,  0.,  0.000018257195229563586,  0.0007682571952295636,  -0.0007682158187528167,   0.0007682571952295636,  -0.0007682158187528167  ]
pyf_kn4 = [  0.,  0.,  0., -0.00001168174301104629,   0.0002883182569889537,  -0.0002883874749225424,   0.0002883182569889537,  -0.0002883874749225424  ]

# beam2.dkp4.final
pxf_dkp4 = [ -0.0010426014142623763,  -0.0010426014142623763,  -0.0010426014142623763,  -0.0010608586094919398,  -0.0003108586094919398,  -0.0017743855955095595,  -0.0003108586094919398,  -0.0017743855955095595  ]
pyf_dkp4 = [  0.00010958182433295981,  0.00010958182433295981,  0.00010958182433295981,  0.00012126356734400609,  0.00042126356734400604, -0.00020203070074449776,  0.00042126356734400604, -0.00020203070074449776 ]

# beam2.dkn4.final
pxf_dkn4 = [  0.0010426014142623763,   0.0010426014142623763,   0.0010426014142623763,   0.0010608586094919398,   0.0018108586094919398,   0.00027438559550955945,  0.0018108586094919398,   0.00027438559550955945 ]
pyf_dkn4 = [ -0.00010958182433295981, -0.00010958182433295981, -0.00010958182433295981, -0.00012126356734400609,  0.00017873643265599388, -0.0003979692992555022,   0.00017873643265599388, -0.0003979692992555022  ]

# -- multipoles

# beam.mp1.final
xf_mp1  = [ 0., 0., 0., 0.0020017506683018874, 0.002039251451642561, -0.0020356721193558517, 0.0020393301102767433, -0.0020357434579643375 ]
pxf_mp1 = [ 0., 0., 0., 0.00007009675756343146, 0.0008216275967063507, -0.0006783111819366863, 0.0008216306774818571, -0.0006783081012716986 ]
yf_mp1  = [ 0., 0., 0., 0.0009956271964705896, 0.0010105357929539218, -0.0010194151082867742, 0.0010105567316816657, -0.0010194541303303404 ]
pyf_mp1 = [ 0., 0., 0., -0.00017508705020997218, 0.00012185301024493958, -0.000477380674232139, 0.0001218468578115431, -0.0004773867853792345 ]
zf_mp1  = [ 0., 0.000049034281550199474, -0.000049037118656134856, 0.00004903383782519922, 0.00004901753480121002, 0.00004901755877446561, -0.00004905393257095815, -0.00004905390847772298 ]

# beam.mn1.final
xf_mn1  = [ 0., 0., 0., 0.0019982493316981127, 0.0020356736483681886, -0.0020392529806509555, 0.0020357449900350244, -0.0020393316423434633 ]
pxf_mn1 = [ 0., 0., 0., -0.00007009675756343146, 0.0006783724032936493, -0.0008216888180633137, 0.0006783693225181429, -0.0008216918987283015 ]
yf_mn1  = [ 0., 0., 0., 0.0010043728035294104, 0.0010194342470475266, -0.0010105549317102746, 0.0010194733084401676, -0.001010575909787066 ]
pyf_mn1 = [ 0., 0., 0., 0.00017508705020997218, 0.0004781469897550604, -0.00012261932576786093, 0.00047815314218845685, -0.0001226132146207654 ]
zf_mn1  = [ 0., 0.000049034281550199474, -0.000049037118656134856, 0.00004903383782519922, 0.00004901754860363581, 0.00004901753120903166, -0.000049053918689922875, -0.00004905393617754518 ]

# beam.mp2.final
xf_mp2  = [ 0., 0., 0., 0.002007971012238037, 0.00223379274561076, -0.002215756602816231, 0.0022342629617779853, -0.0022161863873705577 ]
pxf_mp2 = [ 0., 0., 0., 0.000053193221242483796, 0.0008101764007453585, -0.0006898152587850374, 0.0008101908017093369, -0.0006898008673304005 ]
yf_mp2  = [ 0., 0., 0., 0.000982770378709363, 0.0010708453787873834, -0.0011089004347417377, 0.0010709834375058624, -0.0011091222047735015 ]
pyf_mp2 = [ 0., 0., 0., -0.00011497900515867248, 0.00017277467093106136, -0.0004267286633957622, 0.0001727495363243748, -0.00042675364504917454 ]
zf_mp2  = [ 0., 0.000034265437828406564, -0.00003435659841023283, 0.00003426423650545832, 0.000034165233021345094, 0.000034167350868953886, -0.0000344572059405368, -0.00003445507861381536 ]

# beam.mn2.final
xf_mn2  = [ 0., 0., 0., 0.0019920289877619632, 0.0022157578526597324, -0.0022337939954430387, 0.0022161876382885176, -0.002234264212684649 ]
pxf_mn2 = [ 0., 0., 0., -0.000053193221242483796, 0.0006898235992546416, -0.0008101847412149626, 0.0006898091982906632, -0.0008101991326695995 ]
yf_mn2  = [ 0., 0., 0., 0.0010172296212906371, 0.0011089748602056402, -0.001070919804240822, 0.001109196802203456, -0.001071058034925285 ]
pyf_mn2 = [ 0., 0., 0., 0.00011497900515867248, 0.0004272253290689386, -0.0001732713366042378, 0.00042725046367562513, -0.00017324635495082543 ]
zf_mn2  = [ 0., 0.000034265437828406564, -0.00003435659841023283, 0.00003426423650545832, 0.00003416731826143406, 0.00003416521914531578, -0.0000344551113627096, -0.00003445721987314801 ]

# beam.mp3.final
xf_mp3  = [ 0., 0., 0., 0.002000196794747446, 0.0025622012038446045, -0.0025616754014487235, 0.0025633268746706565, -0.0025627997355042397 ]
pxf_mp3 = [ 0., 0., 0., 5.253107791815821e-7, 0.0007507019247259888, -0.0007492983835403606, 0.0007507023034270174, -0.0007492980053672516 ]
yf_mp3  = [ 0., 0., 0., 0.000999627124171594, 0.0012243000004272947, -0.0012252483969695809, 0.0012247488319113093, -0.001225699561297213 ]
pyf_mp3 = [ 0., 0., 0., -9.953298779579863e-7, 0.0002987312728173831, -0.00030126285858442426, 0.000298730692440487, -0.0003012634347847696 ]
zf_mp3  = [ 0., 1.9534688194108246e-10, -1.9593380261918217e-10, 1.9510986253304015e-10, -2.440585992455014e-7, -2.439488732733258e-7, -2.4542889483776895e-7, -2.453186461569481e-7 ]

# beam.mn3.final
xf_mn3  = [ 0., 0., 0., 0.001999803205252554, 0.0025616752859647196, -0.002562201088360597, 0.0025627996195908705, -0.002563326758757284 ]
pxf_mn3 = [ 0., 0., 0., -5.253107791815821e-7, 0.0007492980752740112, -0.0007507016164596394, 0.0007492976965729827, -0.0007507019946327484 ]
yf_mn3  = [ 0., 0., 0., 0.0010003728758284061, 0.001225250595496266, -0.0012243021989539772, 0.0012257017657931314, -0.001224751036407225 ]
pyf_mn3 = [ 0., 0., 0., 9.953298779579863e-7, 0.00030126872718261685, -0.0002987371414155757, 0.00030126930755951296, -0.0002987365652152303 ]
zf_mn3  = [ 0., 1.9534688194108246e-10, -1.9593380261918217e-10, 1.9510986253304015e-10, -2.4394944850701753e-7, -2.440591687570665e-7, -2.453192240216601e-7, -2.454294669505286e-7 ]

# beam.mp4.final
xf_mp4  = [ 0., 0., 0., 0.0019999767428981424, 0.0036857651724806983, -0.003685864231025714, 0.0036891400174103224, -0.0036892393991604425 ]
pxf_mp4 = [ 0., 0., 0., -2.069365240857758e-8, 0.000749955900651277, -0.0007500440406949199, 0.0007499558451292198, -0.0007500440960746558 ]
yf_mp4  = [ 0., 0., 0., 0.0010000389770261455, 0.0016744006079651998, -0.0016742516246075218, 0.0016757508434368878, -0.001675601396134794 ]
pyf_mp4 = [ 0., 0., 0., 3.468089170817345e-8, 0.00030006647891036176, -0.00029993391686471783, 0.00030006655281196876, -0.0002999338436431643 ]
zf_mp4  = [ 0., 3.140908277834687e-8, -3.1503450227072763e-8, 3.140908186274627e-8, -7.011731101404635e-7, -7.012026794402352e-7, -7.670218208770383e-7, -7.670515526982756e-7 ]

# beam.mn4.final
xf_mn4  = [ 0., 0., 0., 0.0020000232571018577, 0.0036858642969452814, -0.003685765238400266, 0.003689239465372253, -0.0036891400836221327 ]
pxf_mn4 = [ 0., 0., 0., 2.069365240857758e-8, 0.000750044099348723, -0.0007499559593050801, 0.0007500441548707802, -0.0007499559039253443 ]
yf_mn4  = [ 0., 0., 0., 0.0009999610229738545, 0.0016742511798051893, -0.0016744001631628675, 0.0016756009496761398, -0.0016757503969782334 ]
pyf_mn4 = [ 0., 0., 0., -3.468089170817345e-8, 0.0002999335210896382, -0.0003000660831352821, 0.0002999334471880312, -0.00030006615635683566 ]
zf_mn4  = [ 0., 3.140908277834687e-8, -3.1503450227072763e-8, 3.140908186274627e-8, -7.012025955554829e-7, -7.011730261910017e-7, -7.670514683677524e-7, -7.670217364813577e-7 ]

# -- sector bends

# beam.sb1.final
#xf_sb1  = [ 0.,  0.,                      0.,                     4.4834792779340600e-3,  4.5356143287504990e-3, -4.5356143287504990e-3 ]
#pxf_sb1 = [ 0.,  0.,                      0.,                     1.1607821136240948e-1,  1.1776162121669208e-1, -1.1776162121669208e-1 ]
#yf_sb1  = [ 0.,  0.,                      0.,                     1.2400905948673489e-4,  1.3427609296030678e-4, -1.3427609296030678e-4 ]
#pyf_sb1 = [ 0.,  0.,                      0.,                    -2.8691666098954356e-2, -2.8653744321335432e-2,  2.8653744321335432e-2 ]
#zf_sb1  = [ 0.,  4.903428155019947e-5,   -4.903711865613486e-5,  -4.8701323139842656e-5, -5.1605970700562340e-5, -5.1605970700562340e-5 ]


# test individual elements
@testset "ExactTracking" begin
  @testset "Particles" begin

    # ===  D R I F T  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_drift!, (β1, γsq1, 1/βγ1, ld1)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr1 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_drift!, (β2, γsq2, 1/βγ2, ld2)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr2 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_drift!, (β3, γsq3, 1/βγ3, ld3)))
    @test v[:,BeamTracking.XI]  ≈  xf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_dr3 (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] == pxi
    @test v[:,BeamTracking.PYI] == pyi
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_drift!, (β4, γsq4, 1/βγ4, ld4)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1,  gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1, -gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2,  gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2, -gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3,  gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3, -gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4,  gr4 / Bρ4, lq4)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4, -gr4 / Bρ4, lq4)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k1,  kn1 * lk1,  ks1 * lk1, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k1, -kn1 * lk1, -ks1 * lk1, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_dk1 * cos(ra1) / Bρ1
    ks1 = bv_dk1 * sin(ra1) / Bρ1
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk1,  kn1 * lk1,  ks1 * lk1, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk1, -kn1 * lk1, -ks1 * lk1, -1)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k2,  kn2 * lk2,  ks2 * lk2, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k2, -kn2 * lk2, -ks2 * lk2, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn2 = bv_dk2 * cos(ra2) / Bρ2
    ks2 = bv_dk2 * sin(ra2) / Bρ2
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk2,  kn2 * lk2,  ks2 * lk2, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk2, -kn2 * lk2, -ks2 * lk2, -1)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k3,  kn3 * lk3,  ks3 * lk3, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k3, -kn3 * lk3, -ks3 * lk3, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn3 = bv_dk3 * cos(ra3) / Bρ3
    ks3 = bv_dk3 * sin(ra3) / Bρ3
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk3,  kn3 * lk3,  ks3 * lk3, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk3, -kn3 * lk3, -ks3 * lk3, -1)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k4,  kn4 * lk4,  ks4 * lk4, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_k4, -kn4 * lk4, -ks4 * lk4, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn4 = bv_dk4 * cos(ra4) / Bρ4
    ks4 = bv_dk4 * sin(ra4) / Bρ4
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk4,  kn4 * lk4,  ks4 * lk4, -1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.multipole_kick!, (ms_dk4, -kn4 * lk4, -ks4 * lk4, -1)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β1, γsq1, 1/βγ1, 0, ms_m1,  kn1,  ks1, lm1)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β1, γsq1, 1/βγ1, 0, ms_m1, -kn1, -ks1, lm1)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β2, γsq2, 1/βγ2, 0, ms_m2,  kn2,  ks2, lm2)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β2, γsq2, 1/βγ2, 0, ms_m2, -kn2, -ks2, lm2)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β3, γsq3, 1/βγ3, 0, ms_m3,  kn3,  ks3, lm3)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β3, γsq3, 1/βγ3, 0, ms_m3, -kn3, -ks3, lm3)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β4, γsq4, 1/βγ4, 0, ms_m4,  kn4,  ks4, lm4)))
    @test v[:,BeamTracking.XI]  ≈  xf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(IntegrationTracking.dkd_multipole!, (β4, γsq4, 1/βγ4, 0, ms_m4, -kn4, -ks4, lm4)))
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
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_sbend!, (β1, Bρ1, hc1, b_1, ee1, ex1, la1)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb1
    #
    # 1 MeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_sbend!, (β2, Bρ2, hc2, b_2, ee2, ex2, la2)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb2
    #
    # 1 GeV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_sbend!, (β3, Bρ3, hc3, b_3, ee3, ex3, la3)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb3
    #
    # 250 GeV proton
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(Bunch(v).coords, KernelCall(ExactTracking.exact_sbend!, (β4, Bρ4, hc4, b_4, ee4, ex4, la4)))
    @test v[:,BeamTracking.XI]  ≈  xf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb4
    =#

    ###### Exact Sector Bend ##########
    p0c = 10E6
    tilde_m = BeamTracking.massof(Species("electron"))/p0c
    beta_0 = 1/sqrt(1 + tilde_m^2)
    I = SA[1 0 0; 0 1 0; 0 0 1]

    exact_bend_1 =  
      [ 0.8773527130168902E+00  0.4793669072377229E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1225248770305213E+00
       -0.4799049641428075E+00  0.8775825618903755E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.4794255386042034E+00
        0.0000000000000000E+00  0.0000000000000000E+00 0.1000000000000000E+01 0.4999794461108593E+00 0.0000000000000000E+00  0.0000000000000000E+00
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00
       -0.4794255937019158E+00 -0.1222950422117283E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 -0.1925165307287538E-01
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 ]
    g = 1
    L = 0.5
    theta = g * L 
    k0 = 1.001
    test_matrix(exact_bend_1, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))

    exact_bend_2 = 
      [ 0.1139493927324543E+01  0.6225083696592777E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.1231653667943533E-15   
        0.2887206289024770E-15  0.8775825618903689E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.4794255386042027E+00  
        0.0000000000000000E+00  0.0000000000000000E+00 0.1000000000000000E+01 0.5463024898437906E+00 0.0000000000000000E+00  0.0000000000000000E+00  
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00 
       -0.5463024898437888E+00 -0.2984464104095239E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01  0.1302199336067747E-02   
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01 ]
    g = 1
    k0 = 0
    L = 0.5
    theta = g * L
    test_matrix(exact_bend_2, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))

    exact_bend_3 = 
      [ 0.1000000000000000E+01 0.5598925109558526E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1330944687907870E+00  
        0.0000000000000000E+00 0.9999999999999969E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.2784231178942774E-15   
        0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.5186281544969957E+00 0.0000000000000000E+00 0.0000000000000000E+00  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00 0.0000000000000000E+00 
        0.0000000000000000E+00 0.1330944687907864E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.4256655579492583E-01  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01  ]
    k0 = 0.9
    g = 0
    L = 0.5
    theta = g * L
    test_matrix(exact_bend_3, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))

    exact_bend_4 = 
      [ 0.1127528195871212E+01  0.6609770864392946E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.1472939489803990E+00   
        0.1589354646360484E+00  0.9800665778412415E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1986693307950612E+00  
        0.0000000000000000E+00  0.0000000000000000E+00 0.1000000000000000E+01 0.5481504769109217E+00 0.0000000000000000E+00  0.0000000000000000E+00  
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00  
       -0.2474155043455749E+00 -0.2756737519477058E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01  0.7169048328908401E-01   
        0.0000000000000000E+00  0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01  ]
    k0 = -0.8
    g = 0.4
    L = 0.5
    theta = g * L
    test_matrix(exact_bend_4, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))

    exact_bend_5 = 
      [ 0.1283686050523820E+01 0.7804509524043607E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1379043012645095E+00  
        0.2336510053851919E+00 0.9210609940028860E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.3894183423086500E+00   
        0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.5829765210092684E+00 0.0000000000000000E+00  0.0000000000000000E+00   
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00   
        0.5321123724771200E+00 0.4309401889384765E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01  0.8346614392319884E-01  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01   ]
    k0 = 0.6
    g = -0.8
    L = 0.5
    theta = g * L
    test_matrix(exact_bend_5, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))

    exact_bend_6 = 
      [ 0.9005939074669641E+00 0.4858928367401870E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.1117303574601966E+00  
       -0.4314829847437832E+00 0.8775825618903742E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 -0.4794255386042036E+00  
        0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.5022656235734632E+00 0.0000000000000000E+00  0.0000000000000000E+00  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 0.0000000000000000E+00  0.0000000000000000E+00  
        0.4799774672744288E+00 0.1348968216172418E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.1000000000000000E+01 -0.2098595696296651E-01  
        0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00 0.0000000000000000E+00  0.1000000000000000E+01  ]
    k0 = -0.9
    g = -1
    L = 0.5
    theta = g * L
    test_matrix(exact_bend_6, KernelCall(ExactTracking.exact_bend!, (0, 0, theta, g, k0, I, I, tilde_m, beta_0, L)))
  end


  @testset "Utility functions" begin
    dx_rot = -0.1
    dy_rot = -0.1
    dz_rot = 0.2

    sx, cx = sincos(dx_rot/2)
    sy, cy = sincos(dy_rot/2)
    sz, cz = sincos(dz_rot/2)

    W = SA[cx*cy*cz + sx*sy*sz, 
           cy*cz*sx + cx*sy*sz,
           cx*cz*sy - cy*sx*sz,
           -cz*sx*sy + cx*cy*sz]'

    # Test rot_quaternion function
    @test rot_quaternion(dx_rot, dy_rot, dz_rot) ≈ W

    Winv = SA[cx*cy*cz + sx*sy*sz,
              -cy*cz*sx - cx*sy*sz,
              -cx*cz*sy + cy*sx*sz,
              cz*sx*sy - cx*cy*sz]'

    # Test inv_rot_quaternion function
    @test inv_rot_quaternion(dx_rot, dy_rot, dz_rot) ≈ Winv
  end

  @testset "Kernels" begin
    function patch_args(::Type{T}) where {T}
        p0c = T(10e6)
        mc2 = T(BeamTracking.massof(Species("electron")))
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        dt = T(1e-9)
        dx = T(2)
        dy = T(3)
        dz = T(4)
        winv = ExactTracking.inv_rot_quaternion(T(-5),T(6),T(7))
        w31 = 2*(winv[QX]*winv[QZ] - winv[QY]*winv[Q0])
        w32 = 2*(winv[QY]*winv[QZ] + winv[QX]*winv[Q0])
        w33 = 1 - 2*(winv[QX]^2 + winv[QY]^2)
        L = w31*dx + w32*dy + w33*dz
        return beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, winv, L
    end

    function patch_norot_args(::Type{T}) where {T}
        p0c = T(10e6)
        mc2 = T(BeamTracking.massof(Species("electron")))
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
        mc2 = T(BeamTracking.massof(Species("electron")))
        tilde_m = mc2/p0c
        gamsqr_0 = 1 + 1/tilde_m^2
        beta_0 = 1/sqrt(1 + tilde_m^2)
        return beta_0, gamsqr_0, tilde_m, L
    end
    
    function solenoid_args(::Type{T}) where {T}
        L = T(1)
        ks = T(2)
        p0c = T(10e6)
        mc2 = T(BeamTracking.massof(Species("electron")))
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
  end
end
