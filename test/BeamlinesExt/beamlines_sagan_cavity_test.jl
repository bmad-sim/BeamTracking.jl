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

R0 = Beamlines.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out1
  @test b1.coords.q ≈ qout1
  @test b1.t_ref ≈ 6.671281911876738e-9

  ele = sc2
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out2
  @test b1.coords.q ≈ qout2
  @test b1.t_ref ≈ 6.6712819105865366e-9

  ele = sc3
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out3
  @test b1.coords.q ≈ qout3
  @test b1.t_ref == 0

  ele = sc4
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out4
  @test b1.coords.q ≈ qout4
  @test b1.t_ref == 0

  ele = sc5
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out5
  @test b1.coords.q ≈ qout5
  @test b1.t_ref ≈ 6.671281910560153e-9

  ele = sc6
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  if printit; println(b1.coords.v); end
  if printit; println(b1.coords.q); end
  @test b1.coords.v ≈ out6
  @test b1.coords.q ≈ qout6
  @test b1.t_ref ≈ 6.6712819100116676e-9

  ele = sc7
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
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
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  b1.t_ref = t0
  track!(b1, ele)
  @test b1.coords.v ≈ out1
  @test b1.t_ref ≈ t0 + 6.671281911876738e-9

  ele = sc3 # L = 0
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
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
    b = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
    track!(b, ele)
    return b
  end
  @test bs[1].coords.v ≈ bs[2].coords.v
  @test bs[1].coords.q ≈ bs[2].coords.q
  @test bs[1].t_ref ≈ bs[2].t_ref

  # A misaligned cavity (x and y offsets) is the same as shifting the Bunch into the cavity frame
  cavity(sp=species; kw...) = Branch([Marker(E_ref = E0, species_ref = sp), RFCavity(; L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, 
                                     dE_ref = 0.1*E0, tracking_method = SaganCavity(n_cells = 3, L_active = 1.0), Ksol = 0.01, Kn1L = 0.001, kw...)]).beamlines[2].line[1]
  ele = cavity(x_offset = 2e-3, y_offset = -1e-3)
  b1 = Bunch(deepcopy(v1), deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  v2 = deepcopy(v1); v2[1] -= 2e-3; v2[3] += 1e-3
  b2 = Bunch(v2, deepcopy(quat1), species=species, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b2, cavity())
  b2.coords.v[1] += 2e-3; b2.coords.v[3] -= 1e-3
  @test b1.coords.v == b2.coords.v
  @test b1.coords.q == b2.coords.q

  # After a cavity with an energy change, the Bunch reference is the exit reference of the cavity, 
  # so the next Beamline does not ramp the Bunch again (a drift leaves pz unchanged), also for a 
  # species other than the reference species
  for sp in (species, Species("muon"))
    sc = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency = 1e9, dE_ref = 0.1*E0, 
                  tracking_method = SaganCavity(n_cells = 2))
    br = Branch([Beamline([Drift(L = 1.0)]; E_ref = E0, species_ref = species), Beamline([sc]), 
                 Beamline([Drift(L = 1.0)])])
    b = Bunch(deepcopy(v1), deepcopy(quat1), species = sp, 
              p_over_q_ref = br.beamlines[1].p_over_q_ref)
    track!(b, br.beamlines[1]; notify = false)
    track!(b, br.beamlines[2]; notify = false)
    @test b.p_over_q_ref == br.beamlines[2].p_over_q_ref
    pz = b.coords.v[6]
    track!(b, br.beamlines[3]; notify = false)
    @test b.coords.v[6] == pz
  end

  # A Bunch species other than the reference species: a species identical to the electron except 
  # for its name tracks exactly like the electron (with the reference energy change, RF phase and 
  # design clock of a different species)
  e2 = Species("electron2", [getfield(species, n) for n in fieldnames(Species)][2:end]...)
  ele = cavity(Ks0 = 0.001, x_offset = 1e-3, tilt = 0.1)
  bs = map((species, e2)) do sp
    b = Bunch(deepcopy(v1), deepcopy(quat1), species=sp, p_over_q_ref = Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
    track!(b, ele; notify=false)
    return b
  end
  @test bs[1].coords.v == bs[2].coords.v
  @test bs[1].coords.q == bs[2].coords.q
  @test bs[1].t_ref == bs[2].t_ref

  # A proton through a cavity designed for deuterons (same |Bρ|): with dE_ref = 0 and a small voltage 
  # the reference energies are constant, the proton reference particle lags the design deuteron by 
  # s*(1/v_p - 1/v_d) at a kick at s, and gains q*V_k*cos(omega*(t_phi0 + lag)) at each kick. 
  # Afterwards t_ref is the design (deuteron) time and z is referenced to it
  d = Species("deuteron")
  p = Species("proton")
  R = Beamlines.E_to_R(d, 3e9)
  V, f, L, n = 1e3, 5e8, 2.0, 10
  ele = Branch([Marker(p_over_q_ref = R, species_ref = d), RFCavity(L = L, voltage = V, rf_frequency = f, dE_ref = 0.0, 
               tracking_method = SaganCavity(n_cells = n), zero_phase = PhaseRef.AboveTransition, phi0 = 0.3)]).beamlines[2].line[1]
  bp = Bunch(zeros(1, 6); species = p, p_over_q_ref = R)
  track!(bp, ele; notify = false)
  bd = Bunch(zeros(1, 6); species = d, p_over_q_ref = R)
  track!(bd, ele; notify = false)
  @test bp.t_ref == bd.t_ref
  ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
  omega = 2pi*f
  _, L_active = ext.rf_step_calc(n, -1.0, omega, L)
  t_phi0 = ext.rf_phi0_calc(ele.RFParams, d)/omega
  v_p = Beamlines.R_to_v(p, R)
  v_d = Beamlines.R_to_v(d, R)
  P0c = Beamlines.R_to_pc(p, R)
  pz(dE) = sqrt((sqrt(P0c^2 + massof(p)^2) + dE)^2 - massof(p)^2)/P0c - 1
  s_k = [(L - L_active)/2 + k*L_active/n for k in 0:n]
  w_k = [(k == 0 || k == n ? 0.5 : 1.0)/n for k in 0:n]
  dE_slip = sum(V*w_k[k]*cos(omega*(t_phi0 + s_k[k]*(1/v_p - 1/v_d))) for k in 1:n+1)
  dE_no_slip = sum(V*w_k[k]*cos(omega*t_phi0) for k in 1:n+1)
  @test bp.coords.v[6] ≈ pz(dE_slip) rtol = 1e-5
  @test !(bp.coords.v[6] ≈ pz(dE_no_slip))
  @test bp.coords.v[5] ≈ v_p*(L/v_d - L/v_p) rtol = 1e-6

  # Normalized multipoles are normalized by the Beamline p_over_q_ref, so they are the same physical 
  # fields as the unnormalized ones (R times the normalized), also for a positron in an electron cavity
  R_bl = Beamlines.E_to_R(species, E0)
  for sp in (species, Species("positron"))
    vs = map(((; Kn1L = 0.01, Ksol = 0.02), (; Bn1L = 0.01*R_bl, Bsol = 0.02*R_bl), (; Bn1L = -0.01*R_bl, Bsol = -0.02*R_bl))) do kw
      ele = Branch([Marker(E_ref = E0, species_ref = species), RFCavity(; L = 2.0, voltage = 1e3, rf_frequency = 1e9, dE_ref = 0.0, 
                   tracking_method = SaganCavity(n_cells = 3, L_active = 1.0), kw...)]).beamlines[2].line[1]
      b = Bunch(deepcopy(v1), deepcopy(quat1), species=sp, p_over_q_ref = (sp == species ? 1 : -1)*Beamlines.E_to_R(species, ele.E_ref-ele.dE_ref))
      track!(b, ele; notify=false)
      return b.coords.v
    end
    @test vs[1] == vs[2]
    @test maximum(abs.(vs[1] .- vs[3])) > 1e-4
  end

  # The Bunch leaves an accelerating cavity at the exit reference momentum, so a following Drift 
  # does not shift the particles again (on crest with voltage == dE_ref, pz stays ~0)
  let e = Species("electron")
    br = Branch([Marker(E_ref = 1e9, species_ref = e), Drift(L = 1.0),
                 RFCavity(L = 2.0, voltage = 2e8, rf_frequency = 1e9, dE_ref = 2e8, tracking_method = SaganCavity(n_cells = 4)),
                 Drift(L = 1.0)])
    cav_bl = br.beamlines[2]
    b = Bunch([0 0 0 0 0 0.0]; species = e, p_over_q_ref = br.beamlines[1].p_over_q_ref)
    track!(b, br.beamlines[1]; notify = false)
    track!(b, cav_bl.line[1]; notify = false)
    @test b.p_over_q_ref == cav_bl.p_over_q_ref
    pz_cavity = b.coords.v[6]
    @test abs(pz_cavity) < 1e-12
    track!(b, cav_bl.line[2]; notify = false)
    @test b.coords.v[6] == pz_cavity
    @test b.p_over_q_ref == cav_bl.p_over_q_ref
  end
end
