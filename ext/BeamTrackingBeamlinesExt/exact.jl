@inline function pure_patch(tm::Exact, bunch, patchparams, L) 
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  winv = ExactTracking.w_inv_matrix(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(ExactTracking.patch!, (beta_0, gamsqr_0, tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function thick_pure_bsolenoid(tm::Exact, bunch, bm0, L)
  Ks, _ = get_strengths(bm0, L, bunch.Brho_ref)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_solenoid!, (Ks, beta_0, gamsqr_0, tilde_m, L))
end

@inline function drift(tm::Exact, bunch, L)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_drift!, (beta_0, gamsqr_0, tilde_m, L))
end

@inline function thick_bend_pure_bdipole(tm::Exact, bunch, bendparams, bm1, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  w = ExactTracking.w_matrix(0,0,-tilt)
  w_inv = ExactTracking.w_inv_matrix(0,0,-tilt)
  theta = g * L
  Kn0, Ks0 = get_strengths(bm1, L, bunch.Brho_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_bend!, (theta, g, Kn0, w, w_inv, tilde_m, beta_0, L))
end

@inline function thick_pure_bdipole(tm::Exact, bunch, bm1, L)
  Kn0, Ks0 = get_strengths(bm1, L, bunch.Brho_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  I = SA[1 0 0; 0 1 0; 0 0 1]
  return KernelCall(ExactTracking.exact_bend!, (0, 0, Kn0, I, I, tilde_m, beta_0, L))
end

@inline function thick_bend_no_field(tm::Exact, bunch, bendparams, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  w = ExactTracking.w_matrix(0,0,-tilt)
  w_inv = ExactTracking.w_inv_matrix(0,0,-tilt)
  theta = g * L
  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_bend!, (theta, g, 0, w, w_inv, tilde_m, beta_0, L))
end