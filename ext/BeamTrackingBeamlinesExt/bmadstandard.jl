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

@inline function thick_bend_pure_bdipole(tm::BmadStandard, bunch, bendparams, bm1, L)
  g = bendparams.g
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  γ0 = calc_gamma(bunch.species, bunch.Brho_ref)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (g, K0, γ0, G, L))
end

#=@inline function thick_bend_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  g = bendparams.g
  K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
  K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  me2 = bendparams.e2 ≈ 0 ? 0 : K0*tan(bendparams.e2)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (g, k0, γ0, G, L))
end=#

@inline function thick_pure_bquadrupole(tm::BmadStandard, bunch, bm2, L)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  γ0 = calc_gamma(bunch.species, bunch.Brho_ref)
  G = anomalous_moment_of(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_quadrupole!, (K1, γ0, G, L))
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