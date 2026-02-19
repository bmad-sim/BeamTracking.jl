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
out5 = [0.04013970763517789 0.020737423284025818 0.09048039163044863 0.04343792365678535 0.06860289956211306 0.3329611939763494]
out6 = [0.05547643312975879 0.01732154459324357 0.12251750186673921 0.03464880965835598 0.033163420911262356 -0.23655981640919407]
out7 = [0.042245973784445376 0.018306520956017803 0.09275765505991967 0.03600810697674727 0.07171786023122002 0.3544548139812854]

qout1 = [0.999879019609318 -0.013956320891476926 0.006867821712547262 1.6649560539635442e-5]
qout2 = [0.9992430796286998 -0.034793882383567926 0.017396941191783963 0.0]
qout3 = [0.9988275160446326 0.04329982163574608 -0.02164991081787304 0.0]
qout4 = [0.9997215710364805 -0.021105078144390207 0.010552539072195104 0.0]
qout5 = [0.9997167730429309 -0.02057260451935547 0.011964164145470267 2.0401452485922347e-5]
qout6 = [0.9989839935516592 0.040517290217206244 -0.01973139108650261 4.502026884806371e-5]
qout7 = [0.999267804088086 -0.034537890354452774 0.016461727611034916 3.6971580325413343e-5]

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
  @test b1.t_ref ≈  6.925879527784158e-9
end








# dE_ref

