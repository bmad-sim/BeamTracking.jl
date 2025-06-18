@inline function pure_patch(tm::Exact, bunch, patchparams, L) 
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  winv = ExactTracking.w_inv_matrix(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(ExactTracking.patch!, (tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function thick_pure_bsolenoid(tm::Exact, bunch, bm0, L)
  Ks = get_thick_strength(bm0, L, bunch.Brho_ref)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_solenoid!, (Ks, beta_0, gamsqr_0, tilde_m, L))
end

@inline function drift(tm::Exact, bunch, L)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(ExactTracking.exact_drift!, (beta_0, gamsqr_0, tilde_m, L))
end