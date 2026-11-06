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
  cb = (i, coords, cur_s, cur_t_ref, cur_beta_gamma_ref, last_ds_step, last_g, transforms_out!, transforms_in!) -> (i == 1 && push!(ss, cur_s); nothing)
  track_ele(quad, 3, 2; trackkw=(; use_explicit_SIMD=false), callbacks=(cb,))
  ss_sliced = copy(ss)
  empty!(ss)
  track_ele(quad, 1, 6; trackkw=(; use_explicit_SIMD=false), callbacks=(cb,))
  @test ss_sliced ≈ ss
  @test last(ss_sliced) ≈ 2.0

  # Not supported
  @test_throws ErrorException track_ele((; L=2.0), 2, 1; tm=RK4())
  @test_throws ErrorException track_ele((; L=2.0), 2, 1; tm=Implicit())
  @test_throws ErrorException track_ele((; L=2.0, dx=1e-3), 2, 1; tm=Exact())
end
