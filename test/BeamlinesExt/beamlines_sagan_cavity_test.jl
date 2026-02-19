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
quat1 = [1.0  0.0  0.0  0.0]

out1 = [0.999879019609318 -0.013956320891476926 0.006867821712547262 1.6649560539635442e-5]
out2 = [0.9992430796286998 -0.034793882383567926 0.017396941191783963 0.0]
out3 = [0.9988275160446326 0.04329982163574608 -0.02164991081787304 0.0]
out4 = [0.9997215710364805 -0.021105078144390207 0.010552539072195104 0.0]
out5 = [NaN NaN NaN NaN]
out6 = [0.9992383937064061 0.03506860505271674 -0.017112084756381873 4.514408242576583e-5]

R0 = BT.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out1
  println(b1.coords.q)
  @test b1.t_ref ≈ 6.924402749004665e-9

  ele = sc2
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out2
  println(b1.coords.q)
  @test b1.t_ref ≈ 6.8812160436549985e-9

  ele = sc3
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out3
  println(b1.coords.q)
  @test b1.t_ref ≈ 0

  ele = sc4
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out4
  println(b1.coords.q)
  @test b1.t_ref ≈ 0

  ele = sc5
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out5
  println(b1.coords.q)
  @test b1.t_ref ≈ 6.88028440803943e-9

  ele = sc6
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  # println(b1.coords.v)
  @test b1.coords.v ≈ out6
  println(b1.coords.q)
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
