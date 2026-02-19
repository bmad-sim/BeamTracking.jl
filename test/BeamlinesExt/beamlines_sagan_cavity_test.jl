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

include("../lattices/sagan_cavity_lat.jl")
import BeamTracking as BT

v1 = [0.01  0.02  0.03  0.04  0.05  0.1]

out1 = [0.04217873136425864 0.01829391337404973 0.09332500470086044 0.0361983273702167 0.06606835704543228 0.14023826463265046]
out2 = [0.04194241427584856 0.018231932600906685 0.09388482855169711 0.03646386520181337 0.06713366175549244 0.3270142635510056]
out3 = [0.01 0.018382799381485956 0.03 0.03676559876297191 0.04895340854206434 -0.27490000696709455]
out4 = [0.01 0.018048002564646365 0.03 0.03609600512929273 0.05054378625481679 0.20674145046402864]
out5 = [0.04013970763517789 0.020737423284025818 0.09048039163044863 0.04343792365678535 0.06860289956211306 0.3329611939763494]
out6 = [0.05411818209039813 0.017031536764595845 0.11975670804339537 0.03399302609195471 0.03855654036685451 -0.2129528915963015]

R0 = BT.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out1
  @test b1.t_ref ≈ 6.924402749004665e-9

  ele = sc2
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out2
  @test b1.t_ref ≈ 6.8812160436549985e-9

  ele = sc3
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out3
  @test b1.t_ref ≈ 0

  ele = sc4
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out4
  @test b1.t_ref ≈ 0

  ele = sc5
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out5
  @test b1.t_ref ≈ 6.88028440803943e-9

  ele = sc6
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out6
  @test b1.t_ref ≈ 6.862133443570331e-9
end







# Test spin
# Test multipoles
# Test zero length
# standing/traveling waves
# dE_ref
# differing n_steps
# differing n_cell
# non-zero starting orbit
