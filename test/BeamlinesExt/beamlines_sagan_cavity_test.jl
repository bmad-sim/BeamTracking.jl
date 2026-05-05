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

include("../lattices/sagan_cavity_lat.jl")
import BeamTracking as BT
printit = false

v1 = [0.01  0.02  0.03  0.04  0.05  0.1]
quat1 = [1.0  0.0  0.0  0.0]

out1 = [0.04217873136425864 0.01829391337404973 0.09332500470086044 0.0361983273702167 0.06606835704543228 0.14023826463265046]
out2 = [0.04194241427584856 0.018231932600906685 0.09388482855169711 0.03646386520181337 0.06713366175549244 0.3270142635510056]
out3 = [0.01 0.018382799381485956 0.03 0.03676559876297191 0.04895340854206434 -0.27490000696709455]
out4 = [0.01 0.018048002564646365 0.03 0.03609600512929273 0.05054378625481679 0.20674145046402864]
out5 = [0.04024384562264533 0.020905758993788007 0.0906920218736364 0.0437490768778646 0.06860786820748023 0.33296910906572946]
out6 = [0.05547643312975879 0.01732154459324357 0.12251750186673921 0.03464880965835598 0.033163420911262356 -0.23655981640919407]
out7 = [0.039771702242443516 0.020537228663090705 0.08630319003626656 0.040102936717505325 0.07292110767629885 0.35881369057462276]

qout1 = [0.999879019609318 -0.013956320891476926 0.006867821712547262 1.6649560539635442e-5]
qout2 = [0.9992430796286998 -0.034793882383567926 0.017396941191783963 0.0]
qout3 = [0.9988275160446326 0.04329982163574608 -0.02164991081787304 0.0]
qout4 = [0.9997215710364805 -0.021105078144390207 0.010552539072195104 0.0]
qout5 = [0.9997109913949613 -0.020823209473730488 0.012013626268805894 2.0376215140565444e-5]
qout6 = [0.9989839935516592 0.040517290217206244 -0.01973139108650261 4.502026884806371e-5]
qout7 = [0.9998054546005872 -0.01767500106126395 0.008754211736627163 0.00010519157234827785]

R0 = BT.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out1
  @test b1.coords.q ≈ qout1
  @test b1.t_ref ≈ 6.924402749004665e-9

  ele = sc2
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out2
  @test b1.coords.q ≈ qout2
  @test b1.t_ref ≈ 6.8812160436549985e-9

  ele = sc3
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out3
  @test b1.coords.q ≈ qout3
  @test b1.t_ref ≈ 0

  ele = sc4
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out4
  @test b1.coords.q ≈ qout4
  @test b1.t_ref ≈ 0

  ele = sc5
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out5
  @test b1.coords.q ≈ qout5
  @test b1.t_ref ≈ 6.88028440803943e-9

  ele = sc6
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out6
  @test b1.coords.q ≈ qout6
  @test b1.t_ref ≈ 6.862133443570331e-9

  ele = sc7
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out7
  @test b1.coords.q ≈ qout7
  @test b1.t_ref ≈ 6.9251867695848145e-9
end

# dE_ref

