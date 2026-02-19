using Test,
      AtomicAndPhysicalConstants,
      BeamTracking,
      Beamlines,
      JET,
      BenchmarkTools,
      GTPSA,
      StaticArrays,
      ReferenceFrameRotations,
      SIMD

using BeamTracking: Coords, KernelCall, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST, C_LIGHT,
      STATE_LOST_NEG_X, STATE_LOST_POS_X, STATE_LOST_NEG_Y, STATE_LOST_POS_Y, STATE_LOST_PZ, STATE_LOST_Z,
      rot_quaternion, inv_rot_quaternion, atan2, sincu, sinhcu, sincuc, expq, atan2,
      quat_mul, quat_rotate, gaussian_random
using Beamlines: isactive

import BeamTracking as BT

vb1 = [0.1, 0.02, 0.3, 0.04, 0.5, 1.0]
quat1 = [1.0  0.0  0.0  0.0]

out1 = [0.12296022698354489 0.019649290969312984 0.34129839798515305 0.03856589312161067 0.505971927449021 0.7560159043414053]
out2 = [0.12472572480619533 0.019299195161732223 0.3461750920504792 0.037727147222066125 0.5058059437930548 0.7554378514920698]

qout1 = [0.9999707362308135 -0.0038907719842025965 0.0023147865698348053 -0.006166874289958003]
qout2 = [0.9999696325797841 -0.003707386508793316 0.0022676228911355997 -0.0064689326745906595]


m1out = [0.9956867028361168 1.9479825457673512 0.017116161667009182 0.024073951847560882 -2.977620125627317e-5 0.00019628811530733746; -0.004033281596864665 0.9042149609810912 0.005998379632891489 0.01565395420261743 4.32227961792136e-7 -3.306900701841207e-6; -0.004353460968731719 -0.017798273371025013 1.004097771032863 1.9521689562503481 6.19156554670561e-5 -0.0004109015090180842; 0.005999682149132045 -0.0038468072964312344 0.0039678008426686975 0.9121794902618959 1.6613235986401643e-7 -1.1578031529279224e-6; 3.241641652036875e-8 0.00020452630747251006 8.446484776780971e-7 -0.0004073157494565105 0.998715984710241 0.018567966482061833; -6.808669671845504e-8 -1.0307135569267143e-5 -5.106502252988151e-8 2.061427113859003e-5 -0.18723841168108557 0.9059817641041474]

m2out = [0.9690958700629988 1.9454006438640021 0.017517712074800358 0.026066904450224162 -2.637764598683018e-5 0.00020410860260282722; -0.005531100827201287 0.9260219727265299 0.006278551587978387 0.01681122484328427 1.0084655854140437e-5 -1.0663428022519371e-6; -0.004416502802481768 -0.01783115417523026 0.9777493040723166 1.9508899336387555 5.5142669140447805e-5 -0.00042716481805788495; 0.006309160940297908 -0.004150855790860934 0.002860389895862436 0.9345387258963115 -1.9486504708562716e-5 -6.421643263701643e-6; -2.9082509933145476e-6 0.00021282966396674242 6.604384561233256e-6 -0.000423156208271126 0.9986647632724044 0.01854251708208316; 9.922722720281393e-6 4.801412364936898e-6 -2.0131115697102087e-5 -9.671324815522965e-6 -0.187127723723515 0.9060350692247849]

L = 2.0
mass = 1e6
q = -2
E0 = 1e7
dE = 1e6
E1 = E0 + dE
p_over_q = sqrt(E1^2 - mass^2) / (C_LIGHT * q)
t_ref = 0.1 / C_LIGHT
order = SA[0, 1, 2]
Bn = SA[0.01, 0.0001, 0.002] * p_over_q
Bs = SA[0.0, 0.0002, 0.003] * p_over_q
voltage = 4e5
rf_omega = 1e9
t_phi0 = 0.1
a = 0.1

# sagan_cavity_thin!(i, coords::Coords, radiation_damping_on, radiation_fluctuations_on,
#                                      mass, q, E0_ref, dE_ref,
#                                      t_ref, m_order, BnL, BsL, a, q_voltage, rf_omega, t_phi0, L)

# sagan_cavity_thick!(i, coords::Coords, radiation_damping_on, radiation_fluctuations_on, traveling_wave, 
#              mass, q, E0_ref, dE_ref, t_ref, n_cell, m_order, Bn, Bs, 
#              a, q_voltage, rf_omega, t_phi0, L_active, L)



@testset "SaganCavityKernel" begin
  args = (Val{true}(), Val{false}(), mass, q, E0, dE, t_ref, Val{true}(), Val{true}(), order, Bn, Bs, a, q*voltage, rf_omega, t_phi0, L)
  bunch = Bunch(copy(vb1), deepcopy(quat1))
  BT.launch!(bunch.coords, KernelCall(BT.sagan_cavity_zero_L_active!, args))
  @test bunch.coords.v ≈ out1
  @test bunch.coords.q ≈ qout1
  test_matrix(m1out, KernelCall(BT.sagan_cavity_zero_L_active!, args))

  args = (Val{true}(), Val{false}(), Val{false}(), mass, q, E0, dE, t_ref, 2.0, Val{true}(), Val{true}(), order, Bn, Bs, a, q*voltage, rf_omega, t_phi0, L/2, L)
  bunch = Bunch(copy(vb1), deepcopy(quat1))
  BT.launch!(bunch.coords, KernelCall(BT.sagan_cavity_thick!, args))
  @test bunch.coords.v ≈ out2
  @test bunch.coords.q ≈ qout2
  test_matrix(m2out, KernelCall(BT.sagan_cavity_thick!, args))
end
