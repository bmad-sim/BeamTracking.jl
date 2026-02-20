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
printit = false

vb1 = [0.1, 0.02, 0.3, 0.04, 0.5, 1.0]
quat1 = [1.0  0.0  0.0  0.0]


out1 = [0.1 0.018765902124584947 0.3 0.03743180424916989 0.4999561831205719 0.7565003129785745]
out2 = [0.12296022698354489 0.019649290969312984 0.34129839798515305 0.03856589312161067 0.505971927449021 0.7560159043414053]
out3 = [0.1246060556061499 0.019221996565400294 0.3461726582293197 0.037610618298246516 0.5058074766005436 0.7554381718468195]

qout1 = [0.9999945745573493 -0.002928987660974588 0.0015072779663246269 -5.288652091776513e-7]
qout2 = [0.9999707362308135 -0.0038907719842025965 0.0023147865698348053 -0.006166874289958003]
qout3 = [0.9999721220838187 -0.00361182948507389 0.0021972895051256018 -0.006154808021836045]

m1out = [1.0 0.0 0.0 0.0 0.0 0.0; -0.002 0.9082951062292475 0.003 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0; 0.003 0.0 0.002 0.9082951062292475 0.0 0.0; 0.0 0.0 0.0 0.0 1.0004857114354138 0.0; 0.0 0.0 0.0 0.0 -0.18723841643647124 0.9078541510863771]

m2out = [0.9956867028361168 1.9479825457673512 0.017116161667009182 0.024073951847560882 -2.977620125627317e-5 0.00019628811530733746; -0.004033281596864665 0.9042149609810912 0.005998379632891489 0.01565395420261743 4.32227961792136e-7 -3.306900701841207e-6; -0.004353460968731719 -0.017798273371025013 1.004097771032863 1.9521689562503481 6.19156554670561e-5 -0.0004109015090180842; 0.005999682149132045 -0.0038468072964312344 0.0039678008426686975 0.9121794902618959 1.6613235986401643e-7 -1.1578031529279224e-6; 3.241641652036875e-8 0.00020452630747251006 8.446484776780971e-7 -0.0004073157494565105 0.998715984710241 0.018567966482061833; -6.808669671845504e-8 -1.0307135569267143e-5 -5.106502252988151e-8 2.061427113859003e-5 -0.18723841168108557 0.9059817641041474]

m3out = [0.9692196788040458 1.9455452880498731 0.016819827550813207 0.024831899033272847 -2.551516500564454e-5 0.0001989581681828496; -0.005329039569131069 0.926320934913009 0.0059814660860315035 0.015890323703251635 9.864188164944659e-6 -7.802364732510896e-7; -0.004076366857708641 -0.01698866643342548 0.9776420685441772 1.9507752094737267 5.320047183379138e-5 -0.00041541120071264873; 0.006010505804929168 -0.004079846545873315 0.002665322536362557 0.9342640293389503 -1.913789875986876e-5 -6.208161733263679e-6; -2.919895695673319e-6 0.00020693756405682284 6.400402763418292e-6 -0.0004116681416972699 0.9986647648142426 0.01854250744415755; 9.742246939099456e-6 4.651088374924434e-6 -1.9757352824054797e-5 -9.364616457445265e-6 -0.18712772398376687 0.9060350691411453]

L = 2.0
mass = 1e6
q = -2
E0 = 1e7
P0c = sqrt(E0^2 - mass^2)
dE = 1e6
E1 = E0 + dE
p_over_q = sqrt(E1^2 - mass^2) / (C_LIGHT * q)
t_ref = 0.1 / C_LIGHT
order = SA[0, 1, 2]
Bn = SA[0.01, 0.0001, 0.002] * p_over_q
Bs = SA[0.0, 0.0002, 0.003] * p_over_q
voltage = 4e5
L_active = L / 2
gradient = voltage / L_active 
rf_omega = 1e9
t_phi0 = 0.1
a = 0.1

#

@testset "SaganCavityKernel" begin
  args = (Val{true}(), Val{false}(), mass, q, P0c, dE, t_ref, Val{true}(), order, Bn, Bs, a, q*voltage, rf_omega, t_phi0)
  bunch = Bunch(copy(vb1), deepcopy(quat1))
  BT.launch!(bunch.coords, KernelCall(BT.sagan_cavity_zero_L!, args))
  @test bunch.coords.v ≈ out1
  @test bunch.coords.q ≈ qout1
  test_matrix(m1out, KernelCall(BT.sagan_cavity_zero_L!, args), printit = printit)

  args = (Val{true}(), Val{false}(), mass, q, P0c, dE, t_ref, Val{true}(), Val{true}(), order, Bn, Bs, a, q*voltage, rf_omega, t_phi0, L)
  bunch = Bunch(copy(vb1), deepcopy(quat1))
  BT.launch!(bunch.coords, KernelCall(BT.sagan_cavity_zero_L_active!, args))
  @test bunch.coords.v ≈ out2
  @test bunch.coords.q ≈ qout2
  test_matrix(m2out, KernelCall(BT.sagan_cavity_zero_L_active!, args))

  args = (Val{true}(), Val{false}(), Val{false}(), mass, q, P0c, dE, t_ref, 2.0, Val{true}(), Val{true}(), order, Bn, Bs, a, q*gradient, rf_omega, t_phi0, L_active, L)
  bunch = Bunch(copy(vb1), deepcopy(quat1))
  BT.launch!(bunch.coords, KernelCall(BT.sagan_cavity_thick!, args))
  if printit; println(bunch.coords.v); end
  if printit; println(bunch.coords.q); end
  @test bunch.coords.v ≈ out3
  @test bunch.coords.q ≈ qout3
  test_matrix(m3out, KernelCall(BT.sagan_cavity_thick!, args), printit = printit)
end
