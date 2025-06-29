@inline function cavity(tm::BmadStandard, bunch, rfparams, L)
  V = rfparams.voltage
  φ0 = rfparams.phi0
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  if rfparams.harmon_master == true
    wave_number = rfparams.frequency / rfparams.L_ring
  else
    wave_number = rfparams.frequency / (β0 * C_LIGHT)
  end
  p0c = calc_p0c(bunch.species, bunch.Brho_ref)
  return KernelCall(BmadStandardTracking.bmad_cavity!, (V, wave_number, φ0, β0, gamsqr_0, tilde_m, p0c, L))
end

@inline function pure_patch(tm::BmadStandard, bunch, patchparams, L)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  winv = ExactTracking.w_inv_matrix(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(ExactTracking.patch!, (tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function drift(tm::BmadStandard, bunch, L)
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  return KernelCall(BmadStandardTracking.magnus_drift!, (β0, gamsqr_0, tilde_m, L))
end

@inline function thick_pure_bsolenoid(tm::BmadStandard, bunch, bm0, L) 
  Ks = get_thick_strength(bm0, L, bunch.Brho_ref)
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_solenoid!, (Ks, β0, gamsqr_0, tilde_m, G, L))
end

@inline function thick_pure_bdipole(tm::BmadStandard, bunch, bm1, L)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  βγ0 = calc_beta_gamma(bunch.species, bunch.Brho_ref)
  γ0 = sqrt(1 + βγ0^2)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (0, K0, γ0, βγ0, G, L))
end

@inline function thick_bend_pure_bdipole(tm::BmadStandard, bunch, bendparams, bm1, L)
  g = bendparams.g
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  βγ0 = calc_beta_gamma(bunch.species, bunch.Brho_ref)
  γ0 = sqrt(1 + βγ0^2)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (g, K0, γ0, βγ0, G, L))
end

@inline function thick_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  if any(b -> (b > 2 || b == 0), keys(bdict))
    error("BmadStandard does not support thick combined dipoles with higher order multipoles")
  end
  K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
  K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_combined_func!, (0, K0, K1, tilde_m, G, L))
end

@inline function thick_bend_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  if any(b -> (b > 2 || b == 0), keys(bdict))
    error("BmadStandard does not support thick combined dipoles with higher order multipoles")
  end
  g = bendparams.g
  K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
  K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_combined_func!, (g, K0, K1, tilde_m, G, L))
end

@inline function bend_fringe(tm::BmadStandard, bunch, bendparams, bm1, L, upstream)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  e = upstream * bendparams.e1 + (1 - upstream) * bendparams.e2
  return KernelCall(BmadStandardTracking.hwang_edge!, (e, K0, 0, upstream))
end

@inline function bend_fringe(tm::BmadStandard, bunch, bendparams, bm1, bm2, L, upstream)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  e = upstream * bendparams.e1 + (1 - upstream) * bendparams.e2
  return KernelCall(BmadStandardTracking.hwang_edge!, (e, K0, K1, upstream))
end

@inline function thick_pure_bquadrupole(tm::BmadStandard, bunch, bm2, L)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  if abs(K1) < 1e-10
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
    return KernelCall(BmadStandardTracking.magnus_drift!, (β0, gamsqr_0, tilde_m, L))
  else
    βγ0 = calc_beta_gamma(bunch.species, bunch.Brho_ref)
    tilde_m = 1/βγ0
    G = anomalous_moment_of(bunch.species)
    return KernelCall(BmadStandardTracking.magnus_quadrupole!, (K1, βγ0, tilde_m, G, L))
  end
end

@inline function thick_pure_bmultipole(tm::BmadStandard, bunch, bmn, L)
  if bmn.order == 3 # Pure sextupole (No tilt)
    K2 = get_thick_strength(bmn, L, bunch.Brho_ref)
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
    G = anomalous_moment_of(bunch.species)
    return KernelCall(BmadStandardTracking.magnus_thick_sextupole!, (K2, β0, gamsqr_0, tilde_m, G, L))
  
  elseif bmn.order == 4 # Pure octupole (No spin)
    K3N = get_thick_strength(bmn, L, bunch.Brho_ref)*cos(-bmn.order*bmn.tilt)
    K3S = get_thick_strength(bmn, L, bunch.Brho_ref)*sin(-bmn.order*bmn.tilt)
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
    return KernelCall(BmadStandardTracking.thick_octupole!, (K3N, K3S, β0, gamsqr_0, tilde_m, L))

  else
    error("BmadStandard does not support thick multipoles of order greater than octupole")
  end

end