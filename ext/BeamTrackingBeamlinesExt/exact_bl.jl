@inline function fringe_body_fringe(tm::Exact, kc, p_over_q_ref, bunch, ele)
  L = ele.universalparams.L
  isnothing(ele.emfieldparams) || error("Exact is not compatible with elements containing EMFieldParams. Try RK4 instead")
  isnothing(ele.fourpotentialparams) || error("Exact is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.rfparams) || error("Exact is not compatible with elements containing RFParams. Try SolenoidKick or DriftKick instead")
  !isnothing(ele.patchparams) + !isnothing(ele.emultipoleparams) + !isnothing(ele.mapparams) +
    !(isnothing(ele.bmultipoleparams) && isnothing(ele.bendparams)) <= 1 ||
    error("Exact can only handle one of PatchParams, EMultipoleParams, MapParams, or BMultipoleParams/BendParams in a single element")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)
  a = gyromagnetic_anomaly(bunch.species)

  if !isnothing(ele.patchparams)
    # Patch
    params = (beta_0, gamsqr_0, tilde_m, ele.patchparams) # TODO: actual patch kernel params
    return push(kc, integration_launcher(BeamTracking.patch!, params, nothing, nothing, tm, ele.universalparams))
  end

  if !isnothing(ele.mapparams)
    # User-supplied map
    params = (ele.mapparams,) # TODO: actual map kernel params
    return push(kc, integration_launcher(BeamTracking.map!, params, nothing, nothing, tm, ele.universalparams))
  end

  if !isnothing(ele.emultipoleparams)
    # Electric dipole
    L != 0 || error("Exact is not compatible with zero length elements containing EMultipoleParams")
    em = getfield(ele, :emultipoleparams)
    length(em.order) == 1 && first(em.order) == 1 || error("Exact only supports EMultipoleParams with only a dipole component")
    en, es = get_strengths(em, L, p_over_q_ref) # TODO: check electric normalization matches get_strengths
    params = (tilde_m, beta_0, gamsqr_0, a, en[1], es[1]) # TODO: actual electric dipole kernel params
    return push(kc, integration_launcher(BeamTracking.exact_edipole!, params, nothing, nothing, tm, ele.universalparams))
  end

  if isnothing(ele.bmultipoleparams) && isnothing(ele.bendparams)
    # Nothing in the element: marker if zero length, otherwise pure drift
    L != 0 || return kc
    params = (beta_0, gamsqr_0, tilde_m)
    return push(kc, integration_launcher(BeamTracking.exact_drift!, params, nothing, nothing, tm, ele.universalparams))
  end

  if L == 0
    # Thin multipole kick
    isnothing(ele.bendparams) || error("Exact is not compatible with zero length elements containing BendParams")
    mm = ele.bmultipoleparams.order
    first(mm) != 0 || error("Exact is not compatible with zero length elements with a solenoidal multipole component")
    all(ele.bmultipoleparams.integrated) || error("Zero length elements require integrated multipole strengths")
    # get_strengths divides integrated strengths by L, so pass L = 1 to get them back integrated
    kn, ks = get_strengths(ele.bmultipoleparams, one(L), p_over_q_ref)
    params = (tilde_m, beta_0, a, mm, kn, ks) # TODO: actual thin kick kernel params
    return push(kc, integration_launcher(BeamTracking.thin_multipole_kick!, params, nothing, nothing, tm, ele.universalparams))
  end

  q = chargeof(bunch.species)
  mc2 = massof(bunch.species)
  E_ref = mc2/tilde_m/beta_0
  radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing
  fluctuations = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations

  if !isnothing(ele.bendparams)
    g = ele.bendparams.g_ref
    ntilt = -ele.bendparams.tilt_ref
    w = ntilt == 0 ? nothing : rot_quaternion(0, 0, ntilt)
    w_inv = ntilt == 0 ? nothing : inv_rot_quaternion(0, 0, ntilt)
    if isnothing(ele.bmultipoleparams)
      # Bend with no field
      params = (nothing, tilde_m, beta_0, a, g, w, w_inv, 0, SA[0], SA[0], SA[0])
      return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, nothing, nothing, tm, ele.universalparams))
    end
    # Dipole + BendParams
    mm = ele.bmultipoleparams.order
    length(mm) == 1 && first(mm) == 1 || error("Exact only supports BendParams with a pure dipole BMultipoleParams. Try BendKick instead")
    kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref)
    ks[1] == 0 || error("Exact is not currently compatible in an element with both BendParams and a skew dipole multipole")
    Kn0 = kn[1]
    edge_params = (a, tilde_m, beta_0, Kn0, w, w_inv, ele.bendparams.e1, ele.bendparams.e2)
    params = (radiation_params, tilde_m, beta_0, a, g, w, w_inv, Kn0, mm, kn, ks)
    photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, g, w, w_inv, mm, kn, ks) : nothing
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, photon_params, edge_params, tm, ele.universalparams))
  end

  mm = ele.bmultipoleparams.order
  length(mm) == 1 || error("Exact only supports a pure solenoid or pure dipole in a nonzero length element. Try DriftKick, SolenoidKick, MatrixKick or BendKick instead")
  kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref)
  photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, 0, nothing, nothing, mm, kn, ks) : nothing

  if first(mm) == 0
    # Pure solenoid
    edge_params = (a, tilde_m, kn[1], nothing, nothing, nothing, nothing, nothing, nothing)
    params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks)
    return push(kc, integration_launcher(BeamTracking.sks_multipole!, params, photon_params, edge_params, tm, ele.universalparams))
  elseif first(mm) == 1
    # Pure straight dipole, tracked in the dipole frame
    Kn0 = sqrt(kn[1]^2 + ks[1]^2)
    tilt0 = atan(ks[1], kn[1])
    w = tilt0 == 0 ? nothing : rot_quaternion(0, 0, tilt0)
    w_inv = tilt0 == 0 ? nothing : inv_rot_quaternion(0, 0, tilt0)
    edge_params = (a, tilde_m, nothing, kn[1], ks[1], nothing, nothing, nothing, nothing)
    params = (radiation_params, tilde_m, beta_0, a, 0, w, w_inv, Kn0, mm, SA[Kn0], SA[zero(Kn0)])
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, photon_params, edge_params, tm, ele.universalparams))
  else
    error("Exact only supports a pure solenoid or pure dipole in a nonzero length element. Try DriftKick or MatrixKick instead")
  end
end