@inline function cavity(tm::BmadStandard, bunch, rfparams, L)
  V = rfparams.voltage
  φ0 = rfparams.phi0
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  if rfparams.harmon_master == true
    wave_number = rfparams.harmonic_number / rfparams.L_ring
  else
    wave_number = rfparams.rf_frequency / (β0 * C_LIGHT)
  end
  p0c = R_to_pc(bunch.species, bunch.R_ref)
  return KernelCall(BmadStandardTracking.bmad_cavity!, (V, wave_number, φ0, β0, gamsqr_0, tilde_m, p0c, L))
end

@inline function pure_patch(tm::BmadStandard, bunch, patchparams, L)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  winv = ExactTracking.w_inv_quaternion(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(ExactTracking.patch!, (beta_0, gamsqr_0, tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function misalign(tm::BmadStandard, bunch, ap, in)
  tilde_m = massof(bunch.species)/R_to_pc(bunch.species, bunch.R_ref)
  sign = 2*in - 1
  dx = sign * ap.x_offset
  dy = sign * ap.y_offset
  dz = sign * ap.z_offset
  if ap.x_rot != 0 || ap.y_rot != 0 || ap.tilt != 0
    @warn "Rotational misalignments are ignored (currently not supported)"
  end
  winv = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]
  return KernelCall(ExactTracking.patch!, (tilde_m, 0, dx, dy, dz, winv, 0))
end

@inline function drift(tm::BmadStandard, bunch, L)
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  return KernelCall(BmadStandardTracking.magnus_drift!, (β0, gamsqr_0, tilde_m, L))
end

@inline function thick_pure_bsolenoid(tm::BmadStandard, bunch, bm0, L) 
  Ks = get_strengths(bm0, L, bunch.R_ref)
  tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  G = BeamTracking.anom(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_solenoid!, (Ks, β0, gamsqr_0, tilde_m, G, L))
end

@inline function thick_pure_bdipole(tm::BmadStandard, bunch, bm1, L)
  Kn0, Ks0 = get_strengths(bm1, L, bunch.R_ref)
  K0 = sqrt(Kn0^2 + Ks0^2)
  tilt = atan(Ks0, Kn0)
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  βγ0 = R_to_beta_gamma(bunch.species, bunch.R_ref)
  γ0 = sqrt(1 + βγ0^2)
  G = BeamTracking.anom(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (0, 0, 0, K0, w, w_inv, γ0, βγ0, G, L))
end

@inline function thick_bend_pure_bdipole(tm::BmadStandard, bunch, bendparams, bm1, L)
  g = bendparams.g_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  tilt = bendparams.tilt_ref
  Kn0, Ks0 = get_strengths(bm1, L, bunch.R_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot yet be used in a bend")
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  βγ0 = R_to_beta_gamma(bunch.species, bunch.R_ref)
  γ0 = sqrt(1 + βγ0^2)
  G = BeamTracking.anom(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_sbend!, (e1, e2, g, Kn0, w, w_inv, γ0, βγ0, G, L))
end

@inline function thick_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  if any(b -> (b.order > 2 || b.order == 0), bdict)
    error("BmadStandard does not support thick combined dipoles with higher order multipoles")
  end
  kn, ks = get_strengths(bdict, L, bunch.R_ref)
  mm = bdict.order
  k0 = sign(kn[1])*sqrt(kn[1]^2 + ks[1]^2) * (mm[1] == 1)
  k1 = sign(kn[2])*sqrt(kn[2]^2 + ks[2]^2) * (mm[2] == 2)
  if k1 ≈ 0
    thick_pure_bdipole(tm::BmadStandard, bunch, first(bm), L)
  end
  if !(ks[1] ≈ 0 && ks[2] ≈ 0)
    error("BmadStandard does not support tilted combined function magnets")
  end 
  tilde_m = massof(bunch.species)/R_to_pc(bunch.species, bunch.R_ref)
  G = BeamTracking.anom(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_combined_func!, (0, k0, k1, tilde_m, G, L))
end

@inline function thick_bend_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  if any(b -> (b.order > 2 || b.order == 0), bdict)
    error("BmadStandard does not support thick combined dipoles with higher order multipoles")
  end
  g = bendparams.g_ref
  kn, ks = get_strengths(bdict, L, bunch.R_ref)
  mm = bdict.order
  k0 = sign(kn[1])*sqrt(kn[1]^2 + ks[1]^2) * (mm[1] == 1)
  k1 = sign(kn[2])*sqrt(kn[2]^2 + ks[2]^2) * (mm[2] == 2)
  if k1 ≈ 0
    thick_bend_pure_bdipole(tm::BmadStandard, bunch, bendparams, first(bm), L)
  end
  if !(ks[1] ≈ 0 && ks[2] ≈ 0)
    error("BmadStandard does not support tilted combined function magnets")
  end
  tilde_m = massof(bunch.species)/R_to_pc(bunch.species, bunch.R_ref)
  G = BeamTracking.anom(bunch.species)
  return KernelCall(BmadStandardTracking.magnus_combined_func!, (g, k0, k1, tilde_m, G, L))
end

@inline function bend_fringe(tm::BmadStandard, bunch, bendparams, bm1, L, upstream)
  e = upstream * bendparams.e1 + (1 - upstream) * bendparams.e2
  g = bendparams.g_ref
  K0 = get_strengths(bm1, L, bunch.R_ref)
  G = BeamTracking.anom(bunch.species)
  βγ0 = R_to_beta_gamma(bunch.species, bunch.R_ref)
  return KernelCall(BmadStandardTracking.hwang_edge!, (e, g, K0, 0, G, βγ0, upstream))
end

@inline function bend_fringe(tm::BmadStandard, bunch, bendparams, bm1, bm2, L, upstream)
  e = upstream * bendparams.e1 + (1 - upstream) * bendparams.e2
  g = bendparams.g_ref
  K0 = get_strengths(bm1, L, bunch.R_ref)
  K1 = get_strengths(bm2, L, bunch.R_ref)
  G = BeamTracking.anom(bunch.species)
  βγ0 = R_to_beta_gamma(bunch.species, bunch.R_ref)
  return KernelCall(BmadStandardTracking.hwang_edge!, (e, g, K0, K1, G, βγ0, upstream))
end

@inline function thick_pure_bquadrupole(tm::BmadStandard, bunch, bm2, L)
  Kn1, Ks1 = get_strengths(bm2, L, bunch.R_ref)
  K1 = sign(Kn1)*sqrt(Kn1^2 + Ks1^2)
  if abs(K1) < 1e-10
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
    return KernelCall(BmadStandardTracking.magnus_drift!, (β0, gamsqr_0, tilde_m, L))
  else
    βγ0 = R_to_beta_gamma(bunch.species, bunch.R_ref)
    tilde_m = 1/βγ0
    G = BeamTracking.anom(bunch.species)
    return KernelCall(BmadStandardTracking.magnus_quadrupole!, (K1, βγ0, tilde_m, G, L))
  end
end

@inline function thick_pure_bmultipole(tm::BmadStandard, bunch, bmn, L)
  if bmn.order == 3 # Pure sextupole (No tilt)
    Kn2, Ks2 = get_strengths(bmn, L, bunch.R_ref)
    K2 = sign(Kn2)*sqrt(Kn2^2 + Ks2^2)
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
    G = BeamTracking.anom(bunch.species)
    return KernelCall(BmadStandardTracking.magnus_thick_sextupole!, (K2, β0, gamsqr_0, tilde_m, G, L))
  
  elseif bmn.order == 4 # Pure octupole (No spin)
    Kn3, Ks3 = get_strengths(bmn, L, bunch.R_ref)
    tilde_m, gamsqr_0, β0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
    return KernelCall(BmadStandardTracking.thick_octupole!, (Kn3, Ks3, β0, gamsqr_0, tilde_m, L))

  else
    error("BmadStandard does not support thick multipoles of order greater than octupole")
  end

end