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
  R_ref = bunch.R_ref
  mm = bm.order
  knl, ksl = get_integrated_strengths(bm, 0, R_ref)
  params = (mm, knl, ksl, -1)
  if isnothing(bunch.coords.q)
    return KernelCall(ExactTracking.multipole_kick!, params)
  else  
    tilde_m = 1/BeamTracking.R_to_beta_gamma(bunch.species, R_ref)
    return KernelCall(IntegrationTracking.integrate_with_spin_thin!, 
      (ExactTracking.multipole_kick!, params, BeamTracking.anom(bunch.species), 0, tilde_m, mm, knl, ksl))
  end
end

@inline thin_bdipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_pure_bquadrupole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_bquadrupole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_pure_bmultipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)

@inline thin_bmultipole(tm::SplitIntegration, bunch, bm) = thin_pure_bdipole(tm, bunch, bm)


# === Thick elements === #
@inline drift(tm::Union{SplitIntegration,DriftKick}, bunch, L) = drift(Exact(), bunch, L)

@inline function thick_pure_bsolenoid(tm::Union{SplitIntegration,SolenoidKick}, bunch, bm, L) 
  if isnothing(bunch.coords.q)
    return thick_pure_bsolenoid(Exact(), bunch, bm, L)
  else
    R_ref = bunch.R_ref
    tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
    mm = SA[bm.order]
    Ksol, Ksol_skew = get_strengths(bm, L, R_ref)
    kn = SA[Ksol]
    ks = SA[Ksol_skew]
    params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), Ksol, mm, kn, ks)
    return integration_launcher!(IntegrationTracking.sks_multipole!, params, tm, L)
  end
end

@inline function thick_bsolenoid(tm::Union{SplitIntegration,SolenoidKick}, bunch, bm, L) 
  R_ref = bunch.R_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  Ksol = kn[1]
  params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), Ksol, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.sks_multipole!, params, tm, L)
end

@inline function thick_pure_bdipole(tm::DriftKick, bunch, bm, L)
  R_ref = bunch.R_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), mm, kn, ks)
  return integration_launcher!(IntegrationTracking.dkd_multipole!, params, tm, L)
end

@inline thick_bdipole(tm::DriftKick, bunch, bm, L) = thick_pure_bdipole(tm, bunch, bm, L)

@inline function thick_pure_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bm1, L) 
  if isnothing(bunch.coords.q)
    return thick_pure_bdipole(Exact(), bunch, bm1, L)
  else
    R_ref = bunch.R_ref
    tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
    mm = bm1.order
    kn, ks = get_strengths(bm1, L, R_ref)
    k0 = sqrt(kn^2 + ks^2)
    tilt = atan(ks, kn)
    w = ExactTracking.w_quaternion(0,0,tilt)
    w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
    params = (tilde_m, beta_0, BeamTracking.anom(bunch.species), 0, 0, 0, w, w_inv, k0, mm, kn, ks)
    return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
  end
end

@inline function thick_bdipole(tm::BendKick, bunch, bm, L)
  R_ref = bunch.R_ref
  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  k0 = sqrt(kn[1]^2 + ks[1]^2)
  tilt = atan(ks[1], kn[1])
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  params = (tilde_m, beta_0, BeamTracking.anom(bunch.species), 0, 0, 0, w, w_inv, k0, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
end

@inline function thick_bdipole(tm::MatrixKick, bunch, bm, L)
  R_ref = bunch.R_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  k1 = sqrt(kn[2]^2 + ks[2]^2) * (mm[2] == 2)
  if k1 ≈ 0
    return thick_bdipole(DriftKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)
  end
  tilt = (atan(ks[2], kn[2]) / 2) * (mm[2] == 2)
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), w, w_inv, k1, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.mkm_quadrupole!, params, tm, L)
end

@inline function thick_bdipole(tm::SplitIntegration, bunch, bm, L)
  if bm.order[2] == 2
    return thick_bdipole(MatrixKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)
  else
    return thick_bdipole(BendKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)
  end
end

@inline function thick_pure_bquadrupole(tm::Union{SplitIntegration,MatrixKick}, bunch, bm, L)
  R_ref = bunch.R_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  k1 = sqrt(kn^2 + ks^2)
  if k1 ≈ 0
    return thick_pure_bquadrupole(DriftKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)
  end
  tilt = atan(ks, kn) / 2
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), w, w_inv, k1, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.mkm_quadrupole!, params, tm, L)
end

@inline thick_pure_bquadrupole(tm::DriftKick, bunch, bm, L) = 
  thick_pure_bdipole(tm, bunch, bm, L)

@inline function thick_bquadrupole(tm::Union{SplitIntegration,MatrixKick}, bunch, bm, L)
  R_ref = bunch.R_ref
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
  mm = bm.order
  kn, ks = get_strengths(bm, L, R_ref)
  k1 = sqrt(kn[1]^2 + ks[1]^2)
  if k1 ≈ 0
    return thick_bquadrupole(DriftKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)
  end
  tilt = atan(ks[1], kn[1]) / 2
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  params = (beta_0, gamsqr_0, tilde_m, BeamTracking.anom(bunch.species), w, w_inv, k1, mm, kn, ks)
  return integration_launcher!(IntegrationTracking.mkm_quadrupole!, params, tm, L)
end

@inline thick_bquadrupole(tm::DriftKick, bunch, bm, L) = thick_pure_bdipole(tm, bunch, bm, L)

@inline thick_pure_bmultipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L) = 
  thick_pure_bdipole(DriftKick(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bm, L)

@inline thick_bmultipole(tm::Union{SplitIntegration,DriftKick}, bunch, bm, L) = 
  thick_pure_bmultipole(tm, bunch, bm, L)


# =========== BENDING ELEMENTS ============= #
@inline thick_bend_no_field(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, L) = 
  thick_bend_no_field(Exact(), bunch, bendparams, L)

@inline function thick_bend_pure_bdipole(tm::Union{SplitIntegration,BendKick}, bunch, bendparams, bm1, L)
  if isnothing(bunch.coords.q)
    return thick_bend_pure_bdipole(Exact(), bunch, bendparams, bm1, L)
  else
    R_ref = bunch.R_ref
    tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, R_ref)
    g = bendparams.g_ref
    tilt = bendparams.tilt_ref
    e1 = bendparams.e1
    e2 = bendparams.e2
    theta = g * L
    mm = bm1.order
    Kn0, Ks0 = get_strengths(bm1, L, R_ref)
    Ks0 ≈ 0 || error("A skew dipole field cannot yet be used in a bend")
    w = ExactTracking.w_quaternion(0,0,tilt)
    w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
    params = (tilde_m, beta_0, BeamTracking.anom(bunch.species), e1, e2, g, w, w_inv, Kn0, mm, Kn0, Ks0)
    return integration_launcher!(IntegrationTracking.bkb_multipole!, params, tm, L)
  end
end

# =========== PATCH ============= #
@inline pure_patch(tm::SplitIntegration, bunch, patchparams, L)  = 
  pure_patch(Exact(), bunch, patchparams, L)