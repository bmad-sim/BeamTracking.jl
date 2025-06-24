c_light = 299792458.0
degree = π / 180


# ========== define element paameters ==========

# -- drifts
ld1 = 0.15;  # m
ld2 = 0.75;  # m
ld3 = 2.00;  # m
ld4 = 2.00;  # m

# -- quadrupoles
lq1 = 0.05;  # m
lq2 = 0.30;  # m
lq3 = 0.75;  # m
lq4 = 1.20;  # m

gr1 = 0.20;  # T/m
gr2 = 0.60;  # T/m
gr3 = 1.20;  # T/m
gr4 = 3.50;  # T/m

# -- thin-lens kicks
lk1 = 0.05;  # m
lk2 = 0.30;  # m
lk3 = 0.75;  # m
lk4 = 2.25;  # m

ra1 = 15degree
ms_k1  = [ 2, 3, 4, 5, 6 ]
bv_k1  = [ 0.20, 0.18, 0.16, 0.14, 0.12 ];  # T·m^{n-1}
ms_dk1 = [ 1, 2, 3, 4, 5, 6 ]
bv_dk1 = [ 4.e-4, 0.20, 0.18, 0.16, 0.14, 0.12 ];  # T·m^{n-1}

ra2 = 12degree
ms_k2  = [ 2, 3, 4, 5, 6 ]
bv_k2  = [ 0.60, 0.40, 0.30, 0.20, 0.10 ];  # T·m^{n-1}
ms_dk2 = [ 1, 2, 3, 4, 5, 6 ]
bv_dk2 = [ 12.e-4, 0.60, 0.40, 0.30, 0.20, 0.10 ];  # T·m^{n-1}

ra3 = 9degree
ms_k3  = [ 2, 3, 4, 5, 6 ]
bv_k3  = [ 1.20, 1.00, 0.80, 0.60, 0.40 ];  # T·m^{n-1}
ms_dk3 = [ 1, 2, 3, 4, 5, 6 ]
bv_dk3 = [ 0.03, 1.20, 1.00, 0.80, 0.60, 0.40 ];  # T·m^{n-1}

ra4 = 6degree
ms_k4  = [ 2, 3, 4, 5, 6 ]
bv_k4  = [ 3.60, 3.00, 2.40, 1.80, 1.20 ];  # T·m^{n-1}
ms_dk4 = [ 1, 2, 3, 4, 5, 6 ]
bv_dk4 = [ 0.39, 3.60, 3.00, 2.40, 1.80, 1.20 ];  # T·m^{n-1}

# -- multipoles
lm1 = 0.05;  # m
lm2 = 0.30;  # m
lm3 = 0.75;  # m
lm4 = 2.25;  # m
ms_m1 = [ 3, 4, 5, 6 ]
bv_m1 = [ 0.4, 0.3, 0.2, 0.1 ];  # T·m^{n-1}

ra1 = 15degree
ms_m1  = [ 2, 3, 4, 5, 6 ]
bv_m1  = [ 0.20, 0.18, 0.16, 0.14, 0.12 ];  # T·m^{n-1}
ms_dm1 = [ 1, 2, 3, 4, 5, 6 ]
bv_dm1 = [ 4.e-4, 0.20, 0.18, 0.16, 0.14, 0.12 ];  # T·m^{n-1}

ra2 = 12degree
ms_m2  = [ 2, 3, 4, 5, 6 ]
bv_m2  = [ 0.60, 0.40, 0.30, 0.20, 0.10 ];  # T·m^{n-1}
ms_dm2 = [ 1, 2, 3, 4, 5, 6 ]
bv_dm2 = [ 12.e-4, 0.60, 0.40, 0.30, 0.20, 0.10 ];  # T·m^{n-1}

ra3 = 9degree
ms_m3  = [ 2, 3, 4, 5, 6 ]
bv_m3  = [ 1.20, 1.00, 0.80, 0.60, 0.40 ];  # T·m^{n-1}
ms_dm3 = [ 1, 2, 3, 4, 5, 6 ]
bv_dm3 = [ 0.03, 1.20, 1.00, 0.80, 0.60, 0.40 ];  # T·m^{n-1}

ra4 = 6degree
ms_m4  = [ 2, 3, 4, 5, 6 ]
bv_m4  = [ 3.60, 3.00, 2.40, 1.80, 1.20 ];  # T·m^{n-1}
ms_dm4 = [ 1, 2, 3, 4, 5, 6 ]
bv_dm4 = [ 0.39, 3.60, 3.00, 2.40, 1.80, 1.20 ];  # T·m^{n-1}


# ========== define beams ==========

# -- species
#e_minus = Species("electron")
#p_plus =  Species("proton")
#mec2 = massof(e_minus) # 0.51099895069 MeV
#mpc2 = massof(p_plus)  # 938.27208943 MeV
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

# -- initial conditions
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
xf_mp1  = [ 0.,  0.,                      0.,                     2.0038906474613948e-3,  2.0414381618753773e-3, -2.0334840740152908e-3 ]
pxf_mp1 = [ 0.,  0.,                      0.,                     1.5578151070769097e-4,  9.0918337810996800e-4, -5.9070183442735890e-4 ]
yf_mp1  = [ 0.,  0.,                      0.,                     9.9028590464704330e-4,  1.0051011842090445e-3, -1.0248329717192377e-3 ]
pyf_mp1 = [ 0.,  0.,                      0.,                    -3.8895234385991640e-4, -9.5748642016826030e-5, -6.9431180204891080e-4 ]
zf_mp1  = [ 0.,  4.9034281550199474e-5,  -4.9037118656134856e-5,  4.903209153457148e-5,   4.901571516330201e-5, 4.9015774852461346e-5   ]
pzf_mp1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam.mn1.final
xf_mn1  = [ 0.,  0.,                      0.,                     1.9961093525386053e-3,  2.0334869408519317e-3, -2.041441028695591e-3  ]
pxf_mn1 = [ 0.,  0.,                      0.,                    -1.5578151070769097e-4,  5.9081662189003210e-4, -9.092981655726411e-4  ]
yf_mn1  = [ 0.,  0.,                      0.,                     1.0097140953529567e-3,  1.0248688568666960e-3, -1.005137069338170e-3  ]
pyf_mn1 = [ 0.,  0.,                      0.,                     3.8895234385991640e-4,  6.9574864201682590e-4,  9.431180204891087e-5  ]
zf_mn1  = [ 0.,  4.9034281550199474e-5,  -4.9037118656134856e-5,  4.903209153457148e-5,   4.901574824419741e-5,   4.901571596604123e-5  ]
pzf_mn1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]


# test individual elements
@testset "ExactTracking" begin
  #=
  @testset "Particles" begin
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

    # ===  Q U A D R U P O L E  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1,  gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β1, γsq1, 1/βγ1, -gr1 / Bρ1, lq1)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2,  gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β2, γsq2, 1/βγ2, -gr2 / Bρ2, lq2)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3,  gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β3, γsq3, 1/βγ3, -gr3 / Bρ3, lq3)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4,  gr4 / Bρ4, lq4)))
    @test v[:,BeamTracking.XI]  ≈  xf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qf4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qf4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.mkm_quadrupole!, (β4, γsq4, 1/βγ4, -gr4 / Bρ4, lq4)))
    @test v[:,BeamTracking.XI]  ≈  xf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_qd4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_qd4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_qd4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi

    # ===  T H I N - L E N S   K I C K  ===
    #
    # 5 keV electron
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_k1 * cos(ra1) / Bρ1
    ks1 = bv_k1 * sin(ra1) / Bρ1
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k1,  kn1 * lk1,  ks1 * lk1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k1, -kn1 * lk1, -ks1 * lk1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn1 = bv_dk1 * cos(ra1) / Bρ1
    ks1 = bv_dk1 * sin(ra1) / Bρ1
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk1,  kn1 * lk1,  ks1 * lk1)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk1, -kn1 * lk1, -ks1 * lk1)))
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
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k2,  kn2 * lk2,  ks2 * lk2)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k2, -kn2 * lk2, -ks2 * lk2)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn2 = bv_dk2 * cos(ra2) / Bρ2
    ks2 = bv_dk2 * sin(ra2) / Bρ2
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk2,  kn2 * lk2,  ks2 * lk2)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk2, -kn2 * lk2, -ks2 * lk2)))
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
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k3,  kn3 * lk3,  ks3 * lk3)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k3, -kn3 * lk3, -ks3 * lk3)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn3 = bv_dk3 * cos(ra3) / Bρ3
    ks3 = bv_dk3 * sin(ra3) / Bρ3
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk3,  kn3 * lk3,  ks3 * lk3)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk3, -kn3 * lk3, -ks3 * lk3)))
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
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k4,  kn4 * lk4,  ks4 * lk4)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_k4, -kn4 * lk4, -ks4 * lk4)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_kn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    kn4 = bv_dk4 * cos(ra4) / Bρ4
    ks4 = bv_dk4 * sin(ra4) / Bρ4
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk4,  kn4 * lk4,  ks4 * lk4)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2
    v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
    BeamTracking.launch!(BunchView(Bunch(v)), KernelCall(ExactTracking.multipole_kick!, (ms_dk4, -kn4 * lk4, -ks4 * lk4)))
    @test v[:,BeamTracking.XI]  == xi2
    @test v[:,BeamTracking.YI]  == yi2
    @test v[:,BeamTracking.ZI]  == zi2
    @test v[:,BeamTracking.PXI] ≈  pxf_dkn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_dkn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi2


    # ===  M U L T I P O L E  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    kn1 = bv_m1 * cos(ra1) / Bρ1
    ks1 = bv_m1 * sin(ra1) / Bρ1
    #=
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β1, γsq1, 1/βγ1, ms1,  kn1,  ks1, lm1)
    @test v[:,BeamTracking.XI]  ≈  xf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β1, γsq1, 1/βγ1, ms1, -kn1, -ks1, lm1)
    @test v[:,BeamTracking.XI]  ≈  xf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzi
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    kn2 = bv2 * cos(15degree) / Bρ2
    ks2 = bv2 * sin(15degree) / Bρ2
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β2, γsq2, 1/βγ2, ms2, kn2, ks2, lm2)
    @test v[:,BeamTracking.XI]  ≈  xf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mp2
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β2, γsq2, 1/βγ2, ms2, -kn2, -ks2, lm2)
    @test v[:,BeamTracking.XI]  ≈  xf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mn2
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    kn3 = bv3 * cos(15degree) / Bρ3
    ks3 = bv3 * sin(15degree) / Bρ3
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β3, γsq3, 1/βγ3, ms3, kn3, ks3, lm3)
    @test v[:,BeamTracking.XI]  ≈  xf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mp3
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β3, γsq3, 1/βγ3, ms3, -kn3, -ks3, lm3)
    @test v[:,BeamTracking.XI]  ≈  xf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mn3
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    kn4 = bv4 * cos(15degree) / Bρ4
    ks4 = bv4 * sin(15degree) / Bρ4
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β4, γsq4, 1/βγ4, ms4, kn4, ks4, lm4)
    @test v[:,BeamTracking.XI]  ≈  xf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mp4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mp4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mp4
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.dkd_multipole!, v, work, β4, γsq4, 1/βγ4, ms4, -kn4, -ks4, lm4)
    @test v[:,BeamTracking.XI]  ≈  xf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_mn4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_mn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_mn4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_mn4


    # ===  S B E N D  ===
    #
    # 5 keV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.exact_sbend!, v, work, β1, Bρ1, hc1, b_1, ee1, ex1, la1)
    @test v[:,BeamTracking.XI]  ≈  xf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb1  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb1 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb1
    #
    # 1 MeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.exact_sbend!, v, work, β2, Bρ2, hc2, b_2, ee2, ex2, la2)
    @test v[:,BeamTracking.XI]  ≈  xf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb2  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb2 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb2
    #
    # 1 GeV electron
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.exact_sbend!, v, work, β3, Bρ3, hc3, b_3, ee3, ex3, la3)
    @test v[:,BeamTracking.XI]  ≈  xf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb3  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb3 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb3
    #
    # 250 GeV proton
    v = [ xi pxi yi pyi zi pzi ]
    BeamTracking.launch!(ExactTracking.exact_sbend!, v, work, β4, Bρ4, hc4, b_4, ee4, ex4, la4)
    @test v[:,BeamTracking.XI]  ≈  xf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.YI]  ≈  yf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.ZI]  ≈  zf_sb4  (rtol=5.e-13)
    @test v[:,BeamTracking.PXI] ≈  pxf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PYI] ≈  pyf_sb4 (rtol=5.e-13)
    @test v[:,BeamTracking.PZI] == pzf_sb4
    =#
  end
=#
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
        return tilde_m, dt, dx, dy, dz, winv, L
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
        return tilde_m, dt, dx, dy, dz, nothing, L
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
  end
end