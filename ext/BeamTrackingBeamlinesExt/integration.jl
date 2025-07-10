# =========== HELPER FUNCTIONS ============= #
@inline function integration_launcher!(ker, params, tm, L)
  order = tm.order
  ds_step = tm.ds_step
  num_steps = tm.num_steps
  if ds_step < 0
    ds_step = L / num_steps
  else
    num_steps = Int(ceil(L / ds_step))
    ds_step = L / num_steps
  end
  if order == 2
    return KernelCall(IntegrationTracking.order_two_integrator!, (ker, params, ds_step, num_steps, L))
  elseif order == 4
    return KernelCall(IntegrationTracking.order_four_integrator!, (ker, params, ds_step, num_steps, L))
  elseif order == 6
    return KernelCall(IntegrationTracking.order_six_integrator!, (ker, params, ds_step, num_steps, L))
  elseif order == 8
    return KernelCall(IntegrationTracking.order_eight_integrator!, (ker, params, ds_step, num_steps, L))
  end
end


# =========== STRAIGHT ELEMENTS ============= #
# === Thin elements === #
@inline function thin_pure_bdipole(tm::SplitIntegration, bunch, bm)
  brho_0 = bunch.Brho_ref
  mm = bm.order
  knl, ksl = get_integrated_strengths(bm, 0, brho_0)
  return KernelCall(ExactTracking.multipole_kick!, (mm, knl, ksl, 1, 1))
end

@inline thin_bdipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_pure_bquadrupole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_bquadrupole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_pure_bmultipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_bmultipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)


# === Thick elements === #
@inline drift(tm::Union{SplitIntegration,DriftKick}, bunch, L) = drift(Exact(), bunch, L)

@inline thick_pure_bsolenoid(tm::Union{SplitIntegration,SolenoidKick}, bunch, bm, L) = 
  thick_pure_bsolenoid(Exact(), bunch, bm, L)

@inline function thick_bsolenoid(tm::Union{SplitIntegration,SolenoidKick}, bunch, bm, L) 
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  mm = bm.order
  kn, ks = get_strengths(bm, L, brho_0)
  Ksol = sqrt(kn[1]^2 + ks[1]^2)
  params = (beta_0, gamsqr_0, tilde_m, Ksol, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.sks_multipole!, params, tm, L)
end

@inline function thick_pure_bdipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  mm = bm.order
  kn, ks = get_strengths(bm, L, brho_0)
  params = (beta_0, gamsqr_0, tilde_m, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.dkd_multipole!, params, tm, L)
end

@inline thick_bdipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L) = 
  thick_pure_bdipole(tm, bunch, bm, L)

#=
@inline function thick_pure_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bm1, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm1, L, brho_0)
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, 0, b0, 0, 0, L))
end

@inline function thick_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm[1], L, brho_0)
  mm, kn, sn = get_thick_multipoles_no_dipole(bunch, bm, L)
  params = (beta_0, brho_0, 0, b0, 0, 0, mm, kn, sn)
  return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
end
=#
@inline function thick_pure_bquadrupole(tm::Union{SplitIntegration,MatrixKick}, bunch, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  mm = bm.order
  kn, ks = get_strengths(bm, L, brho_0)
  k1 = sqrt(kn[1]^2 + ks[1]^2)
  tilt = bm.tilt[1]
  w = ExactTracking.w_matrix(0,0,-tilt)
  w_inv = ExactTracking.w_inv_matrix(0,0,-tilt)
  params = (beta_0, gamsqr_0, tilde_m, w, w_inv, k1, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.mkm_quadrupole!, params, tm, L)
end

@inline thick_pure_bquadrupole(tm::DriftKick, bunch, bm, L) = 
  thick_pure_bdipole(tm, bunch, bm, L)

@inline thick_bquadrupole(tm::Union{SplitIntegration,MatrixKick}, bunch, bm, L) = 
  thick_pure_bquadrupole(tm, bunch, bm, L)

@inline thick_bquadrupole(tm::DriftKick, bunch, bm, L) = thick_pure_bdipole(tm, bunch, bm, L)

@inline thick_pure_bmultipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L) = 
  thick_pure_bdipole(DriftKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)

@inline thick_bmultipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L) = 
  thick_pure_bmultipole(tm, bunch, bm, L)


# =========== BENDING ELEMENTS ============= #
#=@inline function thick_bend_no_field(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, hc, 0, e1, e2, L))
end

@inline function thick_bend_pure_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm1, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm1, L, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  return KernelCall(ExactTracking.exact_sbend!, (beta_0, brho_0, hc, b0, e1, e2, L))
end

@inline function thick_bend_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  b0 = brho_0 * get_thick_strength(bm.bdict[1], L, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  mm, kn, sn = get_thick_multipoles_no_dipole(bunch, bm, L)
  params = (beta_0, brho_0, hc, b0, e1, e2, mm, kn, sn)
  return integration_launcher!(EIntegrationTracking.bkb_multipole!, params, tm, L)
end

@inline function thick_bend_pure_bquadrupole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm2, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  K1 = get_thick_strength(bm2, L, brho_0)
  mm = @SArray [bm2.order]
  kn = @SArray [K1 * cos(-bm2.order*bm2.tilt)]
  sn = @SArray [K1 * sin(-bm2.order*bm2.tilt)]
  params = (beta_0, brho_0, hc, 0, e1, e2, mm, kn, sn)
  return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
end

@inline function thick_bend_bquadrupole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm, L)
  brho_0 = bunch.Brho_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, brho_0)
  hc, e1, e2 = bendparams.g, bendparams.e1, bendparams.e2
  mm, kn, sn = get_thick_multipoles(bunch, bm, L)
  params = (beta_0, brho_0, hc, 0, e1, e2, mm, kn, sn)
  return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
end

@inline thick_bend_pure_bmultipole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm1, L) = thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm1, L)

@inline thick_bend_bmultipole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bdict, L) = thick_bend_bquadrupole(tm, bunch, bendparams, bdict, L)
=#

# =========== PATCH ============= #
@inline pure_patch(tm::SplitIntegration, bunch, patchparams, L)  = 
  pure_patch(Exact(), bunch, patchparams, L)