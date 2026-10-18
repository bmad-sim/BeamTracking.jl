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


@testset "SaganCavityUniversal" begin
  # The reference time is accumulated through a SaganCavity (not reset)
  t0 = 1e-6
  ele = sc1
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  b1.t_ref = t0
  track!(b1, ele)
  @test b1.coords.v ≈ out1
  @test b1.t_ref ≈ t0 + 6.671281911876738e-9

  ele = sc3 # L = 0
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  b1.t_ref = t0
  track!(b1, ele)
  @test b1.coords.v ≈ out3
  @test b1.t_ref ≈ t0

  # Errors for parameter groups not supported by SaganCavity name the element
  sc_bend = RFCavity(name = "sc_bend", L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, g_ref = 0.01,
                     tracking_method = SaganCavity(n_cells = 2))
  sc_patch = RFCavity(name = "sc_patch", L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, dx = 1e-3,
                      tracking_method = SaganCavity(n_cells = 2))
  for (ele, pg) in ((sc_bend, "BendParams"), (sc_patch, "PatchParams"))
    bl = Beamline([ele], species_ref = species, E_ref = E0)
    b1 = Bunch(deepcopy(v1), deepcopy(quat1), species = species)
    @test_throws "SaganCavity Tracking through element $(ele.name) with $pg is undefined" track!(b1, bl)
  end

  # With the BendParams in ignore_params, the cavity is tracked as if it has no BendParams
  sc_bend = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, dE_ref = 0.1*E0, g_ref = 0.01,
                     tracking_method = SaganCavity(n_cells = 2), ignore_params = [BendParams])
  sc_nobend = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, dE_ref = 0.1*E0,
                       tracking_method = SaganCavity(n_cells = 2))
  bs = map((sc_bend, sc_nobend)) do sc
    ele = Branch([Marker(E_ref = E0, species_ref = species), sc]).beamlines[2].line[1]
    b = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
    track!(b, ele)
    return b
  end
  @test bs[1].coords.v ≈ bs[2].coords.v
  @test bs[1].coords.q ≈ bs[2].coords.q
  @test bs[1].t_ref ≈ bs[2].t_ref
end

@testset "SaganCavity reference momentum" begin
  # Issue 367: after the cavity the bunch must be referenced to the cavity exit energy so the
  # next element does not shift pz a second time.
  for nc in (4, 0)
    cav = RFCavity(L = 1.0, voltage = 1e8, rf_frequency = 1.3e9, dE_ref = 1e8, phi0 = 0.0,
            zero_phase = PhaseRef.Accelerating, tracking_method = SaganCavity(n_cells = nc), traveling_wave = true)
    br = Branch([Marker(E_ref = 1e8, species_ref = species), cav, Drift(L = 0.5), Drift(L = 0.5)])
    eles = collect(Iterators.flatten(bl.line for bl in br.beamlines))
    b = Bunch(zeros(1, 6), species = species, p_over_q_ref = BT.E_to_R(species, 1e8))
    for ele in eles
      track!(b, ele)
      @test abs(b.coords.v[1,6]) < 1e-14
    end
    @test b.p_over_q_ref ≈ eles[2].p_over_q_ref
  end

  # A bunch arriving with a reference different from E0_ref = E_ref - dE_ref is shifted to E0_ref
  # first: the result must match tracking the same particle given relative to E0_ref.
  sc = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, dE_ref = 0.1*E0, phi0 = 0.1,
                tracking_method = SaganCavity(n_cells = 3))
  ele = Branch([Marker(E_ref = E0, species_ref = species), sc]).beamlines[2].line[1]
  E_in = ele.E_ref - ele.dE_ref
  R_in = BT.E_to_R(species, E_in)
  R_b = BT.E_to_R(species, 1.05*E_in)
  v_b = copy(v1)
  v_b[1,2] *= R_in/R_b;  v_b[1,4] *= R_in/R_b;  v_b[1,6] = (1 + v1[1,6])*R_in/R_b - 1
  b_ref = Bunch(copy(v1), species = species, p_over_q_ref = R_in)
  b_mis = Bunch(v_b, species = species, p_over_q_ref = R_b)
  track!(b_ref, ele)
  track!(b_mis, ele)
  @test b_mis.coords.v ≈ b_ref.coords.v
  @test b_mis.p_over_q_ref ≈ ele.p_over_q_ref
  @test b_mis.t_ref ≈ b_ref.t_ref

  # The particle time used to evaluate time dependent parameters at the entrance must use the
  # entrance reference E0_ref (RefState beta_gamma_enter): a voltage that depends on time must give
  # the same result as a constant voltage equal to its value at the particle's entrance time.
  # Low energy so that the particle velocity differs strongly between E0_ref and E_ref.
  V(t) = 3e6*(1 + 1e9*t)
  mk(volt) = Branch([Marker(E_ref = 1e6, species_ref = species),
                     RFCavity(L = 1.0, voltage = volt, rf_frequency = 1e9, dE_ref = 2e6, phi0 = 0.1,
                              tracking_method = SaganCavity(n_cells = 3))]).beamlines[2].line[1]
  ele_t = mk(3e6*(1 + 1e9*Time()))
  R_in = BT.E_to_R(species, 1e6)
  t0 = 2e-9
  K = (1 + v1[1,6])*BT.R_to_beta_gamma(species, R_in)
  t_particle = t0 - v1[1,5]/BT.beta_gamma_to_v(K)
  ele_c = mk(V(t_particle))
  b_t = Bunch(copy(v1), species = species, p_over_q_ref = R_in, t_ref = t0)
  b_c = Bunch(copy(v1), species = species, p_over_q_ref = R_in, t_ref = t0)
  track!(b_t, ele_t)
  track!(b_c, ele_c)
  @test b_t.coords.v ≈ b_c.coords.v
  @test b_t.t_ref ≈ b_c.t_ref

  # Ramped (time dependent) reference energy: equivalent to a constant reference equal to the
  # ramped value at the bunch entrance time.
  E_ramp = TimeDependentParam(t -> E0*(1 + 1e6*t), false)
  t0 = 1e-6
  ele_r = Branch([Marker(E_ref = E_ramp, species_ref = species), deepcopy(sc)]).beamlines[2].line[1]
  ele_c = Branch([Marker(E_ref = E0*(1 + 1e6*t0), species_ref = species), deepcopy(sc)]).beamlines[2].line[1]
  @test ele_r.p_over_q_ref(t0) ≈ ele_c.p_over_q_ref
  R_in = BT.E_to_R(species, ele_c.E_ref - ele_c.dE_ref)
  b_r = Bunch(copy(v1), species = species, p_over_q_ref = R_in, t_ref = t0)
  b_c = Bunch(copy(v1), species = species, p_over_q_ref = R_in, t_ref = t0)
  track!(b_r, ele_r)
  track!(b_c, ele_c)
  @test b_r.coords.v ≈ b_c.coords.v
  @test b_r.p_over_q_ref ≈ ele_c.p_over_q_ref
  @test b_r.t_ref ≈ b_c.t_ref
end
