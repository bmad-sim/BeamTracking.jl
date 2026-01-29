@inline function pure_patch(tm::Exact, bunch, patchparams, L) 
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  winv = inv_rot_quaternion(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(BeamTracking.patch!, (beta_0, gamsqr_0, tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function thick_pure_bsolenoid(tm::Exact, bunch, bm0, L)
  Ksol, _ = get_strengths(bm0, L, bunch.p_over_q_ref)
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_solenoid!, (Ksol, beta_0, gamsqr_0, tilde_m, L))
end

@inline function drift(tm::Exact, bunch, L)
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_drift!, (beta_0, gamsqr_0, tilde_m, L))
end

@inline function thick_bend_pure_bdipole(tm::Exact, bunch, bendparams, bm1, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  Kn0, Ks0 = get_strengths(bm1, L, bunch.p_over_q_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_bend_with_rotation!, (e1, e2, theta, 0, g, Kn0, w, w_inv, tilde_m, beta_0, L))
end

@inline function thick_pure_bdipole(tm::Exact, bunch, bm1, L)
  Kn0, Ks0 = get_strengths(bm1, L, bunch.p_over_q_ref)
  Kn = sqrt(Kn0^2 + Ks0^2)
  tilt = atan2(Ks0, Kn0)
  w = rot_quaternion(0,0,tilt)
  w_inv = inv_rot_quaternion(0,0,tilt)
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_bend_with_rotation!, (0, 0, 0, 0, 0, Kn, w, w_inv, tilde_m, beta_0, L))
end

@inline function thick_bend_no_field(tm::Exact, bunch, bendparams, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_curved_drift!, (e1, e2, theta, g, w, w_inv, gyromagnetic_anomaly(bunch.species), tilde_m, beta_0, L))
end

@inline pure_map(tm::Exact, bunch, mapparams, L) = KernelCall(BeamTracking.map!, (mapparams.transport_map, L))

#=
@inline function bend_entrance_fringe(tm::Exact, bunch, bendparams, bmp, L)
  e1 = bendparams.e1
  tilt = bendparams.tilt_ref
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  mm = bmp.order
  kn, ks = get_strengths(bmp, L, bunch.p_over_q_ref)
  (mm[1] == 1 && ks[1] ≈ 0) || error("A skew dipole field cannot be used in an exact bend")
  Kn0 = ifelse(mm[1] == 1, kn[1], 0)
  return KernelCall(BeamTracking.linear_bend_fringe!, (e1, Kn0, w, w_inv))
end

@inline function bend_exit_fringe(tm::Exact, bunch, bendparams, bmp, L)
  e2 = bendparams.e2
  tilt = bendparams.tilt_ref
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  mm = bmp.order
  kn, ks = get_strengths(bmp, L, bunch.p_over_q_ref)
  (mm[1] == 1 && ks[1] ≈ 0) || error("A skew dipole field cannot be used in an exact bend")
  Kn0 = ifelse(mm[1] == 1, kn[1], 0)
  return KernelCall(BeamTracking.linear_bend_fringe!, (e2, Kn0, w, w_inv))
end
=#
