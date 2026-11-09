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

# A slice of a SaganCavity tracks the part of the one cavity between two cuts (equal length slices): 
# the kicks (set by n_cells), fringes and drifts stay where they are, a drift that crosses a cut is 
# split, and its half kicks stay at its ends. So slicing never changes the cavity: n_slices > 1 tracks
# like n_slices = 1 up to round-off, and bitwise if every cut is on a kick.
@testset "SaganCavitySlices" begin
  ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
  e = Species("electron")
  E = 1e10
  vs0 = [1e-3 1e-4 -1e-3 2e-4 1e-3 1e-4; -2e-3 0 1e-3 -1e-4 0 -2e-4]
  qs0 = [1.0 0 0 0; 1.0 0 0 0]
  sliced_cavity(n_slices; kw...) = Branch([Marker(E_ref = E, species_ref = e), RFCavity(; L = 2.0, voltage = 0.2*E, 
                    rf_frequency = 1e9, dE_ref = 0.1*E, n_slices, kw...)]).beamlines[2].line[1]
  function track_sliced(n_slices, sp; callbacks = (), trackkw = (;), kw...)
    c = sliced_cavity(n_slices; kw...)
    R = (sp == e ? 1 : -1)*Beamlines.E_to_R(e, c.E_ref - c.dE_ref)
    b = Bunch(copy(vs0), copy(qs0); species = sp, p_over_q_ref = R, callbacks)
    track!(b, c; notify = false, trackkw...)
    return b
  end
  slice_cases = [
    (; tracking_method = SaganCavity(n_cells = 4)),
    (; tracking_method = SaganCavity(n_cells = 3, L_active = 1.0), Ksol = 0.01, Kn1L = 0.001, x_offset = 1e-3, y_rot = 1e-3, tilt = 0.1),
    (; tracking_method = SaganCavity(n_cells = 0, L_active = 1.2), traveling_wave = false, phi0 = 0.1),
    (; tracking_method = SaganCavity(n_cells = 3, L_active = 0.0), zero_phase = PhaseRef.AboveTransition, Kn2L = 0.1),
    (; tracking_method = SaganCavity(), Ks0 = 0.001, traveling_wave = true),
    (; tracking_method = SaganCavity(n_cells = 5, L_active = 1.9), fringe_at = Fringe.EntranceEnd),
  ]
  for kw in slice_cases, sp in (e, Species("positron"))
    b1 = track_sliced(1, sp; kw...)
    for n_slices in (2, 3, 4, 7)
      b = track_sliced(n_slices, sp; kw...)
      @test maximum(abs.(b.coords.v .- b1.coords.v)) < 1e-14
      @test maximum(abs.(b.coords.q .- b1.coords.q)) < 1e-14
      @test b.t_ref ≈ b1.t_ref rtol = 1e-15
      @test b.p_over_q_ref == b1.p_over_q_ref
      @test all(b.coords.state .== b1.coords.state)
    end
  end

  # Kicks at L_active/4 steps: every cut is on a kick for n_slices = 2 and 4
  for n_slices in (2, 4), sp in (e, Species("positron"))
    b1 = track_sliced(1, sp; slice_cases[1]...)
    b = track_sliced(n_slices, sp; slice_cases[1]...)
    @test b.coords.v == b1.coords.v
    @test b.coords.q == b1.coords.q
    @test b.t_ref == b1.t_ref
  end

  # Negative controls: one more cell, or no fringes, changes the result far more than the tolerance
  b1 = track_sliced(1, e; slice_cases[2]...)
  @test maximum(abs.(track_sliced(3, e; slice_cases[2]..., tracking_method = SaganCavity(n_cells = 4, L_active = 1.0)).coords.v .- b1.coords.v)) > 1e-8
  @test maximum(abs.(track_sliced(3, e; slice_cases[2]..., fringe_at = Fringe.NoEnd).coords.v .- b1.coords.v)) > 1e-8
  # fringe_at is honored, also when the fringes are inside a slice
  for fa in (Fringe.NoEnd, Fringe.EntranceEnd, Fringe.ExitEnd)
    @test maximum(abs.(track_sliced(3, e; slice_cases[2]..., fringe_at = fa).coords.v .- track_sliced(1, e; slice_cases[2]..., fringe_at = fa).coords.v)) < 1e-14
  end
  @test track_sliced(1, e; slice_cases[2]..., fringe_at = Fringe.EntranceEnd).coords.v != b1.coords.v
  @test track_sliced(1, e; slice_cases[2]..., fringe_at = Fringe.ExitEnd).coords.v != b1.coords.v

  # The windows of the slices partition the cavity: each kick and fringe once, the drift lengths add 
  # up to L, the reference energy at a slice exit is that at the next slice entrance, and the 
  # reference times add up to that of the whole cavity
  for (n_cells, L_active) in ((4, 2.0), (3, 1.0), (0, 1.2), (0, 0.0), (7, 1.9), (1, 2.0)), n_slices in (1, 2, 3, 5, 8)
    L = 2.0
    s_cut(k) = k == n_slices ? L : k*(L/n_slices)
    ws = [ext.sagan_window(n_cells, L_active, L, (k-1)*(L/n_slices), s_cut(k), k == 1, k == n_slices, true, true) for k in 1:n_slices]
    n_kicks = (L_active == 0 || n_cells == 0) ? 1 : n_cells + 1
    @test ws[1].ka == 0
    @test ws[end].kb == n_kicks - 1
    @test all(ws[k+1].ka == ws[k].kb + 1 for k in 1:n_slices-1)
    @test count(w -> w.fringe_in, ws) == (L_active == 0 ? 0 : 1)
    @test count(w -> w.fringe_out, ws) == (L_active == 0 ? 0 : 1)
    L_outer = L_active == 0 ? L/2 : (L - L_active)/2
    L_cell = n_cells == 0 ? L_active/2 : L_active/n_cells
    drift_len(on, full, len, nominal) = on ? (full ? nominal : len) : 0.0
    total = sum(ws) do w
      drift_len(w.outer_in, w.outer_in_full, w.outer_in_len, L_outer) + drift_len(w.outer_out, w.outer_out_full, w.outer_out_len, L_outer) +
      drift_len(w.pre, w.pre_full, w.pre_len, L_cell) + drift_len(w.post, w.post_full, w.post_len, L_cell) + 
      (L_active == 0 || n_cells == 0 ? 0.0 : max(w.kb - w.ka, 0)*L_cell)
    end
    @test total ≈ L rtol = 1e-14
    # Each half kick of a drift is in exactly one slice
    n_end_kicks = sum(w -> w.outer_in_kick_in + w.outer_in_kick_out + w.outer_out_kick_in + w.outer_out_kick_out + 
                           w.pre_kick_in + w.pre_kick_out + w.post_kick_in + w.post_kick_out + 2*max(w.kb - w.ka, 0)*(n_cells > 0), ws)
    @test n_end_kicks == 4 + (L_active == 0 ? 0 : (n_cells == 0 ? 4 : 2*n_cells))
    refs = [ext.sagan_window_ref(e, E, 0.1*E, n_cells, L_active, L, w) for w in ws]
    @test refs[1][2] == E - 0.1*E
    @test refs[end][3] == E
    @test all(refs[k+1][2] == refs[k][3] for k in 1:n_slices-1)
    w1 = ext.sagan_window(n_cells, L_active, L, 0.0, L, true, true, true, true)
    @test sum(r -> r[1], refs) ≈ ext.sagan_window_ref(e, E, 0.1*E, n_cells, L_active, L, w1)[1] rtol = 1e-14
  end

  # Callbacks: one per slice and particle (n_slices*n_steps, see beamlines_slices_test.jl)
  n_cb = Ref(0)
  cb = ((i, coords, s, t, bg, ds, g, o!, in!) -> (n_cb[] += 1; nothing),)
  track_sliced(3, e; callbacks = cb, trackkw = (; use_explicit_SIMD = false), slice_cases[1]...)
  @test n_cb[] == 3*size(vs0, 1)
end
