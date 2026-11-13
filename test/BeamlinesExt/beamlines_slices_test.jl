using Random

@testset "Slices" begin
  Random.seed!(0)
  e = Species("electron")
  R0 = -10.0
  v0 = 1e-3 .* randn(50, 6)

  # Each slice of an element with n_slices = k is tracked like an element of length L/k, with the
  # fringes only at the element ends, so the element must track exactly like k such elements
  # (including the ramping and the collective effects recomputed at the start of each slice)
  function test_slices_vs_elements(kw, k; tm=MatrixKick(), ramp=false, species=e, q=nothing, weight=nothing, v=v0, trackkw=(;))
    sliced = LatElement(; tracking_method=tm, n_slices=k, n_steps=3, kw...)
    vkw = haskey(kw, :voltage) ? (; voltage=kw.voltage/k) : (;)
    pieces = [LatElement(; tracking_method=tm, n_steps=3, kw..., L=kw.L/k, vkw...,
                         fringe_at=(j == 1 ? Fringe.EntranceEnd : (j == k ? Fringe.ExitEnd : Fringe.NoEnd))) for j in 1:k]
    R = ramp ? R0*(1 + 1e4*Time()) : R0
    b1 = Bunch(copy(v), deepcopy(q), weight; p_over_q_ref=(species == e ? R0 : -R0), species=species)
    b2 = Bunch(copy(v), deepcopy(q), weight; p_over_q_ref=(species == e ? R0 : -R0), species=species)
    track!(b1, Beamline([sliced]; p_over_q_ref=R, species_ref=e); notify=false, trackkw...)
    track!(b2, Beamline(pieces; p_over_q_ref=R, species_ref=e); notify=false, trackkw...)
    @test b1.coords.v == b2.coords.v
    @test b1.coords.q == b2.coords.q
    @test b1.t_ref == b2.t_ref
    @test b1.p_over_q_ref == b2.p_over_q_ref
  end

  rf = (; L=4.0, voltage=3e6, rf_frequency=5e8, zero_phase=PhaseRef.AboveTransition)
  test_slices_vs_elements((; L=2.0, Kn1=0.3, Ks2=0.5), 3)
  test_slices_vs_elements((; L=2.0, Kn1=0.3, Ks2=0.5), 3; ramp=true)
  test_slices_vs_elements((; L=2.0, Kn1=0.3, Ks2=0.5), 3; ramp=true, trackkw=(; ramp_update_each_particle=true))
  test_slices_vs_elements((; L=2.0, g=0.1, Kn0=0.1, e1=0.05, e2=0.03), 4; tm=BendKick())
  test_slices_vs_elements((; L=1.5, Ksol=0.3, Kn1=0.1), 3; tm=SolenoidKick(), q=repeat([1.0 0 0 0], 50))
  test_slices_vs_elements((; L=0.8, Kn2=1.3), 2; tm=DriftKick())
  test_slices_vs_elements((; L=3.0), 3; tm=Exact())
  test_slices_vs_elements(rf, 2; tm=DriftKick())
  test_slices_vs_elements(rf, 2; tm=DriftKick(), ramp=true)
  test_slices_vs_elements(rf, 2; tm=DriftKick(), species=Species("positron")) # RF slip
  test_slices_vs_elements((; L=20.0, Kn1=0.01, ibs_model=IBSModel.GaussianDampingOnly), 4; weight=1e10, v=1e-5 .* randn(200, 6))

  # A misaligned element with 2 slices is the same body as 2 elements of length L/2 with the same 
  # rotation, each with its center offset by the rotated half-length: the center of each half (at 
  # c = -L/4 and L/4 in the body frame) is at offset + W*(0, 0, c), with W the rotation of the 
  # misalignment. With x/y offsets and tilt this is bitwise; with z_offset, x_rot and y_rot the 
  # offsets of the halves and the lever arms in the kernels (lengths ~ L) are rounded differently, 
  # so they agree to a few ulps of L
  function test_misaligned_slices_vs_halves(al, test)
    L = 20.0
    kw = (; Kn1=0.01, ibs_model=IBSModel.GaussianDampingOnly)
    wz = BeamTracking.quat_rotate((0.0, 0.0, 1.0), BeamTracking.rot_quaternion(al.x_rot, al.y_rot, al.tilt))
    sliced = LatElement(; L=L, tracking_method=MatrixKick(), n_steps=3, n_slices=2, kw..., al...)
    halves = [LatElement(; L=L/2, tracking_method=MatrixKick(), n_steps=3, kw..., fringe_at=fa, x_rot=al.x_rot, y_rot=al.y_rot, tilt=al.tilt,
                         x_offset=al.x_offset + c*wz[1], y_offset=al.y_offset + c*wz[2], z_offset=al.z_offset + c*(wz[3] - 1)) 
              for (c, fa) in ((-L/4, Fringe.EntranceEnd), (L/4, Fringe.ExitEnd))]
    v = 1e-5 .* randn(200, 6)
    b1 = Bunch(copy(v), nothing, 1e10; p_over_q_ref=R0, species=e)
    b2 = Bunch(copy(v), nothing, 1e10; p_over_q_ref=R0, species=e)
    track!(b1, Beamline([sliced]; p_over_q_ref=R0, species_ref=e); notify=false)
    track!(b2, Beamline(halves; p_over_q_ref=R0, species_ref=e); notify=false)
    @test test(b1.coords.v, b2.coords.v)
    @test b1.t_ref == b2.t_ref
  end
  test_misaligned_slices_vs_halves((; x_offset=1e-3, y_offset=-2e-3, z_offset=0.0, x_rot=0.0, y_rot=0.0, tilt=0.3), ==)
  test_misaligned_slices_vs_halves((; x_offset=1e-3, y_offset=-2e-3, z_offset=3e-3, x_rot=2e-3, y_rot=-1e-3, tilt=0.3), (a, b) -> maximum(abs.(a .- b)) <= 4*eps(20.0))

  # Same for a misaligned bend. The misalignment is the rigid motion x -> C + d + W*(x - C) in the 
  # untilted frame at the arc midpoint, with C the chord midpoint, d the offsets and W the rotation. 
  # The untilted frame at the arc midpoint of each half (sigma = -L/4, L/4 from the arc midpoint) 
  # is y -> P + Q*y (reference orbit P, rotation Q by g*sigma about (sin(tilt_ref), -cos(tilt_ref), 0)), 
  # so each half has the rotation Q'*W*Q and the offsets Q'*(C + d + W*(P - C) - P) - C_h + Q'*W*Q*C_h
  qmat(q) = hcat(collect(BeamTracking.quat_rotate((1.0, 0, 0), q)), collect(BeamTracking.quat_rotate((0, 1.0, 0), q)), collect(BeamTracking.quat_rotate((0, 0, 1.0), q)))
  chord_mid(L, g, st, ct) = (f = -L^2*g*BeamTracking.one_cos_norm(L*g/2)/4; [f*ct, f*st, 0.0])
  function misaligned_bend_halves_test(tilt_ref)
    L, g = 4.0, 0.1
    al = (; x_offset=1e-3, y_offset=-2e-3, z_offset=3e-3, x_rot=2e-3, y_rot=-1e-3, tilt=0.3)
    st, ct = sincos(tilt_ref)
    W = qmat(BeamTracking.rot_quaternion(al.x_rot, al.y_rot, al.tilt))
    C = chord_mid(L, g, st, ct)
    halves_al = map((-L/4, L/4)) do sigma
      P = [-ct, -st, 0.0]*(1 - cos(g*sigma))/g + [0, 0, 1.0]*sin(g*sigma)/g
      Q = qmat((cos(g*sigma/2), sin(g*sigma/2)*st, -sin(g*sigma/2)*ct, 0.0))
      Ch = chord_mid(L/2, g, st, ct)
      Wh = Q'*W*Q
      dh = Q'*(C + [al.x_offset, al.y_offset, al.z_offset] + W*(P - C) - P) - Ch + Wh*Ch
      # Angles of Wh = R_y(y_rot)*R_x(x_rot)*R_z(tilt) (see rot_quaternion)
      x_rot = asin(-Wh[2,3])
      y_rot = atan(Wh[1,3], Wh[3,3])
      M = qmat(BeamTracking.rot_quaternion(x_rot, 0, 0))'*qmat(BeamTracking.rot_quaternion(0, y_rot, 0))'*Wh
      (; x_offset=dh[1], y_offset=dh[2], z_offset=dh[3], x_rot, y_rot, tilt=atan(M[2,1], M[1,1]))
    end
    kw = (; g=g, Kn0=g, tilt_ref=tilt_ref, e1=0.05, e2=0.03, Kn2=0.5, ibs_model=IBSModel.GaussianDampingOnly)
    sliced = LatElement(; L=L, tracking_method=BendKick(), n_steps=3, n_slices=2, kw..., al...)
    halves = [LatElement(; L=L/2, tracking_method=BendKick(), n_steps=3, kw..., halves_al[j]..., fringe_at=(j == 1 ? Fringe.EntranceEnd : Fringe.ExitEnd)) for j in 1:2]
    v = 1e-5 .* randn(200, 6)
    b1 = Bunch(copy(v), nothing, 1e10; p_over_q_ref=R0, species=e)
    b2 = Bunch(copy(v), nothing, 1e10; p_over_q_ref=R0, species=e)
    track!(b1, Beamline([sliced]; p_over_q_ref=R0, species_ref=e); notify=false)
    track!(b2, Beamline(halves; p_over_q_ref=R0, species_ref=e); notify=false)
    @test maximum(abs.(b1.coords.v .- b2.coords.v)) <= 4*eps(L)
    @test b1.t_ref == b2.t_ref
  end
  misaligned_bend_halves_test(0.0)
  misaligned_bend_halves_test(0.4)

  # Without ramping or collective effects, n_slices = k is the same as n_slices = 1 with k times
  # the steps: fringes are only at the element ends, apertures are only checked at the element
  # ends, and misalignments transform out and back in at each slice end
  function track_ele(kw, k, n_steps; tm=MatrixKick(), ramp=false, trackkw=(;), callbacks=())
    ele = LatElement(; tracking_method=tm, n_slices=k, n_steps=n_steps, kw...)
    b = Bunch(copy(v0); p_over_q_ref=R0, species=e, callbacks=callbacks)
    track!(b, Beamline([ele]; p_over_q_ref=(ramp ? R0*(1 + 1e4*Time()) : R0), species_ref=e); notify=false, trackkw...)
    return b
  end
  quad = (; L=2.0, Kn1=0.3, Ks2=0.5)
  for fringe_at in (Fringe.BothEnds, Fringe.EntranceEnd, Fringe.ExitEnd, Fringe.NoEnd)
    @test track_ele((; quad..., fringe_at), 3, 3).coords.v == track_ele((; quad..., fringe_at), 1, 9).coords.v
  end
  @test track_ele(quad, 3, 3).t_ref ≈ track_ele(quad, 1, 9).t_ref
  misaligned = (; x_offset=1e-3, y_offset=-1e-3, z_offset=1e-3, x_rot=1e-3, y_rot=2e-3, tilt=0.1)
  @test track_ele((; quad..., misaligned...), 3, 3).coords.v ≈ track_ele((; quad..., misaligned...), 1, 9).coords.v atol=1e-15
  bend = (; L=2.0, g=0.1, Kn0=0.1, e1=0.05, e2=0.03, misaligned...)
  @test track_ele(bend, 4, 3; tm=BendKick()).coords.v ≈ track_ele(bend, 1, 12; tm=BendKick()).coords.v atol=1e-15
  aperture = (; x1_limit=-1.5e-3, x2_limit=1.5e-3, y1_limit=-1.5e-3, y2_limit=1.5e-3, aperture_at=ApertureAt.BothEnds)
  b1 = track_ele((; quad..., aperture...), 3, 3)
  b2 = track_ele((; quad..., aperture...), 1, 9)
  @test b1.coords.state == b2.coords.state
  @test any(!=(STATE_ALIVE), b1.coords.state) && any(==(STATE_ALIVE), b1.coords.state)
  @test b1.coords.v == b2.coords.v

  # The reference energy is ramped at the start of each slice
  b1 = track_ele(quad, 3, 3; ramp=true)
  b2 = track_ele(quad, 1, 9; ramp=true)
  @test !(b1.coords.v ≈ b2.coords.v)
  @test b1.p_over_q_ref ≈ R0*(1 + 1e4*(2/3*2/Beamlines.R_to_v(e, R0))) rtol=1e-6 # p_over_q_ref at the last slice
  @test b2.p_over_q_ref == R0
  b1 = track_ele(quad, 3, 3; ramp=true, trackkw=(; ramp_update_each_particle=true))
  b2 = track_ele(quad, 1, 9; ramp=true, trackkw=(; ramp_update_each_particle=true))
  @test b1.p_over_q_ref == b2.p_over_q_ref
  @test !(b1.coords.v ≈ b2.coords.v)

  # Callbacks get s from the element entrance
  ss = Float64[]
  cb = (i, coords, cur_s, cur_t_ref, cur_step, cur_beta_gamma_ref, last_ds_step, last_g, transforms_out!, transforms_in!) -> (i == 1 && push!(ss, cur_s); nothing)
  track_ele(quad, 3, 2; trackkw=(; use_explicit_SIMD=false), callbacks=(cb,))
  ss_sliced = copy(ss)
  empty!(ss)
  track_ele(quad, 1, 6; trackkw=(; use_explicit_SIMD=false), callbacks=(cb,))
  @test ss_sliced ≈ ss
  @test last(ss_sliced) ≈ 2.0

  # RK4 and Implicit evaluate the fields at s (and the element-local time, which is computed from 
  # s) from the element entrance, so an element with s- and t-dependent fields tracks the same 
  # with n_slices = k as with n_slices = 1 and k times the steps, and the same as k elements with 
  # the fields shifted by s0 and s0/(beta_0*c). If s and t restarted in each slice it would not
  rk_field(x, y, s, t, p) = (p[1]*x*sin(s) + p[2]*cos(1e8*t), p[1]*y*cos(s), 0.0*x, p[2]*y*sin(s), p[2]*x*cos(s) + p[1]*sin(1e8*t), 0.0*x)
  function implicit_potential(x, y, s, t, p)
    a0 = p
    potential = (a0*BeamTracking.C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), a0*cosh(x)*sin(y)*cos(s)*sinh(t), a0*sinh(x)*cosh(y)*sin(s)*cos(t), a0*cos(x)*sinh(y)*cosh(s)*sin(t))
    jac = (a0*BeamTracking.C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -a0*BeamTracking.C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), a0*BeamTracking.C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), a0*BeamTracking.C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
           a0*cos(s)*sin(y)*sinh(t)*sinh(x),  a0*cos(s)*cos(y)*cosh(x)*sinh(t), -a0*cosh(x)*sin(s)*sin(y)*sinh(t),  a0*cos(s)*cosh(t)*cosh(x)*sin(y),
           a0*cos(t)*cosh(x)*cosh(y)*sin(s),  a0*cos(t)*sin(s)*sinh(x)*sinh(y),  a0*cos(s)*cos(t)*cosh(y)*sinh(x), -a0*cosh(y)*sin(s)*sin(t)*sinh(x),
          -a0*cosh(s)*sin(t)*sin(x)*sinh(y), a0*cos(x)*cosh(s)*cosh(y)*sin(t), a0*cos(x)*sin(t)*sinh(s)*sinh(y), a0*cos(t)*cos(x)*cosh(s)*sinh(y))
    return potential, jac
  end
  bg = Beamlines.R_to_beta_gamma(e, R0)
  beta_0 = bg/sqrt(1 + bg^2)
  shifted(f, s0) = (x, y, s, t, p) -> f(x, y, s + s0, t + s0/(beta_0*BeamTracking.C_LIGHT), p)
  rk4_ele(f=rk_field; L=1.2, kw...) = LatElement(; L=L, Kn1=0.2, em_field=f, em_field_params=(0.01, 0.02), em_field_normalized=true, tracking_method=RK4(), kw...)
  implicit_ele(f=implicit_potential; L=1.2, kw...) = LatElement(; L=L, four_potential=f, four_potential_params=1e-3, four_potential_normalized=true, tracking_method=Implicit(), kw...)
  v8 = 1e-3 .* randn(8, 6)
  q8 = repeat([1.0 0 0 0], 8)
  function track_eles(eles; q=nothing, callbacks=(), trackkw=(;))
    b = Bunch(copy(v8), deepcopy(q); p_over_q_ref=R0, species=e, callbacks=callbacks)
    track!(b, Beamline(collect(eles); p_over_q_ref=R0, species_ref=e); notify=false, trackkw...)
    return b
  end
  maxdiff(a, b) = maximum(abs.(a .- b))
  for (mk, f, kw) in ((rk4_ele, rk_field, (;)), (implicit_ele, implicit_potential, (;)), (implicit_ele, implicit_potential, (; rad_damping=true)))
    qq = mk === implicit_ele ? q8 : nothing
    sliced = track_eles([mk(; n_slices=3, n_steps=3, kw...)]; q=qq)
    unsliced = track_eles([mk(; n_steps=9, kw...)]; q=qq)
    L3 = 1.2/3 # Length of each slice
    pieces = track_eles([mk(shifted(f, (j - 1)*L3); L=L3, n_steps=3, kw..., fringe_at=(j == 1 ? Fringe.EntranceEnd : (j == 3 ? Fringe.ExitEnd : Fringe.NoEnd))) for j in 1:3]; q=qq)
    restarting = track_eles([mk(; L=L3, n_steps=3, kw...) for _ in 1:3]; q=qq)
    @test maxdiff(sliced.coords.v, unsliced.coords.v) <= 1e-16
    @test maxdiff(sliced.coords.v, pieces.coords.v) <= 1e-16
    @test sliced.t_ref == pieces.t_ref
    @test maxdiff(restarting.coords.v, unsliced.coords.v) > 1e-5
    isnothing(qq) || @test maxdiff(sliced.coords.q, unsliced.coords.q) <= 1e-16
  end

  # Inside the element the scalar potential is in general nonzero at the slice ends (Implicit and 
  # the Exact electric dipole convert with it there, as for callbacks), so callbacks in a misaligned 
  # element see the same global coordinates with and without slices
  function global_callback_coords(ele)
    rec = Vector{Float64}[]
    cb = (i, coords, cur_s, cur_t_ref, step, bg, ds, g, transforms_out!, transforms_in!) -> begin
      if i == 1
        transforms_out!(i, coords)
        push!(rec, [cur_s, coords.v[i,1], coords.v[i,2], coords.v[i,5], coords.v[i,6]])
        transforms_in!(i, coords)
      end
      nothing
    end
    track_eles([ele]; callbacks=(cb,), trackkw=(; use_explicit_SIMD=false))
    return reduce(hcat, rec)
  end
  @test global_callback_coords(implicit_ele(; n_slices=3, n_steps=2, misaligned...)) ≈ global_callback_coords(implicit_ele(; n_steps=6, misaligned...)) atol=1e-15
  @test global_callback_coords(rk4_ele(; n_slices=3, n_steps=2, misaligned...)) ≈ global_callback_coords(rk4_ele(; n_steps=6, misaligned...)) atol=1e-15

  # Exact electric dipole: the potential is -kE*x in the frame of the field, so a tilted dipole with 
  # the particle rotated with it is the same as the untilted dipole, also for the momentum at the 
  # slice end and for callbacks (the pz below is the same scalar in both)
  rotated(v, a) = (c = cos(a); sn = sin(a); [c*v[1]-sn*v[3]  c*v[2]-sn*v[4]  sn*v[1]+c*v[3]  sn*v[2]+c*v[4]  v[5]  v[6]])
  function elsep_midpoint_pz(kw, v0, k, n)
    pz = Float64[]
    cb = (i, coords, cur_s, cur_t_ref, step, bg, ds, g, transforms_out!, transforms_in!) -> begin
      if abs(cur_s - 1.0) < 1e-12
        transforms_out!(i, coords)
        push!(pz, coords.v[i,6])
        transforms_in!(i, coords)
      end
      nothing
    end
    ele = LatElement(; L=2.0, tracking_method=Exact(), n_slices=k, n_steps=n, fringe_at=Fringe.NoEnd, kw...)
    b = Bunch(copy(v0); p_over_q_ref=R0, species=e, callbacks=(cb,))
    track!(b, Beamline([ele]; p_over_q_ref=R0, species_ref=e); notify=false, use_explicit_SIMD=false)
    return b.coords.v, pz
  end
  v1 = [3e-3 1e-4 -2e-3 2e-4 1e-4 1e-4]
  etilt = 0.7
  v_untilted, pz_untilted = elsep_midpoint_pz((; En0=2e6), v1, 1, 2)
  for (k, n) in ((1, 2), (2, 1))
    v_tilted, pz_tilted = elsep_midpoint_pz((; En0=2e6*cos(etilt), Es0=2e6*sin(etilt)), rotated(v1, -etilt), k, n)
    @test rotated(v_tilted, etilt) ≈ v_untilted atol=1e-17
    @test pz_tilted ≈ pz_untilted rtol=1e-14
  end

  # Not supported
  @test_throws ErrorException track_ele((; L=2.0, voltage=1e6, rf_frequency=5e8), 2, 1; tm=SaganCavity())
  @test_throws ErrorException track_ele((; L=2.0, dx=1e-3), 2, 1; tm=Exact())
  # Elements tracked in one step error for n_steps > 1, so n_slices*n_steps is always the number of
  # steps (and of callbacks, which e.g. SciBmad's twiss relies on)
  count_cb(n) = ((i, coords, s, t, step, bg, ds, g, o!, in!) -> (n[] += 1; nothing),)
  for kw in ((; L=2.0, dx=1e-3), (; L=2.0, transport_map=(v, q) -> (v, q)))
    @test_throws "n_steps > 1" track_ele(kw, 1, 2; tm=Exact())
    ncb = Ref(0)
    track_ele(kw, 1, 1; tm=Exact(), callbacks=count_cb(ncb), trackkw=(; use_explicit_SIMD=false)) # one call per particle
    @test ncb[] == size(v0, 1)
  end
  sagan(n_steps) = Beamline([Drift(L=0.5), RFCavity(L=2.0, voltage=1e6, rf_frequency=5e8, tracking_method=SaganCavity(n_cells=3),
                             n_steps=n_steps)]; p_over_q_ref=R0, species_ref=e).line[2]
  @test_throws "n_steps > 1" track!(Bunch(copy(v0); p_over_q_ref=R0, species=e), sagan(2); notify=false)
  ncb = Ref(0)
  track!(Bunch(copy(v0); p_over_q_ref=R0, species=e, callbacks=count_cb(ncb)), sagan(1); notify=false, use_explicit_SIMD=false)
  @test ncb[] == size(v0, 1)
end
