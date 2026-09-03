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

out1 = [0.04400975583062183 0.019078953928035665 0.09749725528525972 0.038135350406844216 0.048566293508714706 0.10438207893767194]
out2 = [0.044131463345117555 0.01833333333017374 0.09826292669023509 0.03666666666034748 0.04853795083994871 0.15126623838903783]
out3 = [0.01 0.018461538459062357 0.03 0.036923076918124714 0.049999999987821946 -0.11790579365526294]
out4 = [0.01 0.018181818177698334 0.03 0.03636363635539667 0.050000000009851144 0.10606758069187741]
out5 = [0.04335160449531273 0.019509688087400362 0.09770093134315469 0.040881719689982984 0.04857563116670182 0.15383255215875477]
out6 = [0.04978157815779977 0.018558663846092498 0.11025218958580006 0.03708466852871902 0.04797989591851569 -0.04606025892622582]
out7 = [0.0436921489049126 0.019570972953271013 0.09469762610273985 0.038335723907434234 0.04865287895447131 0.1584594807959742]

qout1 = [0.9969082114110158 0.07015209076904136 -0.03539231914414494 0.0002931341345360551]
qout2 = [0.9882219020700481 0.13687212212589134 -0.06843606106294567 0.0]
qout3 = [0.9841416507357511 -0.15865739513079358 0.07932869756539679 0.0]
qout4 = [0.9943358106063602 0.09506332940011167 -0.047531664700055834 0.0]
qout5 = [0.9957322482368662 0.07718724537569223 -0.05056522513866194 0.0016052972134234728]
qout6 = [0.9907268064380208 -0.12175040740909787 0.060305653940766195 0.0006792669694700435]
qout7 = [0.9945900636739683 0.09519284014707527 -0.040806404330128124 -0.00798534854536703]

R0 = BT.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out1
  @test b1.coords.q ≈ qout1
  @test b1.t_ref ≈ 6.671281911876738e-9

  ele = sc2
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out2
  @test b1.coords.q ≈ qout2
  @test b1.t_ref ≈ 6.6712819105865366e-9

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
  @test b1.t_ref ≈ 6.671281910560153e-9

  ele = sc6
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out6
  @test b1.coords.q ≈ qout6
  @test b1.t_ref ≈ 6.6712819100116676e-9

  ele = sc7
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out7
  @test b1.coords.q ≈ qout7
  @test b1.t_ref ≈ 6.6712819118982366e-9
end

# dE_ref

