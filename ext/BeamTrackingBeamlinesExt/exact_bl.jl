@inline function fringe_body_fringe(tm::Exact, kc, bunch, ele)
  L = ele.universalparams.L
  p_over_q_ref = ele.ref.p_over_q_ref
  isnothing(ele.emfieldparams) || error("Exact is not compatible with elements containing EMFieldParams. Try RK4 instead")
  isnothing(ele.fourpotentialparams) || error("Exact is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.rfparams) || error("Exact is not compatible with elements containing RFParams. Try SolenoidKick or DriftKick instead")
  !isnothing(ele.patchparams) + !isnothing(ele.emultipoleparams) + !isnothing(ele.mapparams) +
    !(isnothing(ele.bmultipoleparams) && isnothing(ele.bendparams)) <= 1 ||
    error("Exact is only compatible with one of PatchParams, EMultipoleParams, MapParams, or BMultipoleParams/BendParams")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)
  a = gyromagnetic_anomaly(bunch.species)

  if !isnothing(ele.patchparams)
    # Patch
    ele.universalparams.n_steps == 1 || error("Elements with n_steps > 1 are not compatible with PatchParams (a patch is one step). Try n_steps = 1 instead")
    pp = ele.patchparams
    winv = inv_rot_quaternion(pp.dx_rot, pp.dy_rot, pp.dz_rot)
    return push(kc, make_kernel_call(BeamTracking.patch!, (beta_0, gamsqr_0, tilde_m, pp.dt, pp.dx, pp.dy, pp.dz, winv, L)))
  end

  if !isnothing(ele.mapparams)
    # User-supplied map, which is in the coordinates of the reference species
    ele.universalparams.n_steps == 1 || error("Elements with n_steps > 1 are not compatible with MapParams (a map is one step). Try n_steps = 1 instead")
    bunch.species == ele.ref.species_ref || error("Exact is not compatible with elements containing MapParams when the Bunch species ($(bunch.species)) differs from the reference species ($(ele.ref.species_ref))")
    return push(kc, make_kernel_call(BeamTracking.map!, (ele.mapparams.transport_map, ele.mapparams.transport_map_params, L)))
  end

  if !isnothing(ele.emultipoleparams)
    # Electric dipole
    L != 0 || error("Exact is not compatible with zero length elements containing EMultipoleParams")
    em = ele.emultipoleparams
    length(em.order) == 1 && first(em.order) == 1 || error("Exact only supports EMultipoleParams with only a dipole component")
    en, es = get_e_strengths(em, L)
    kn = en[1]
    ks = es[1]
    kE = sqrt(kn^2 + ks^2)/(p_over_q_ref * C_LIGHT)
    tilt = atan(ks, kn)
    w = tilt == 0 ? nothing : rot_quaternion(0, 0, tilt)
    w_inv = tilt == 0 ? nothing : inv_rot_quaternion(0, 0, tilt)
    q = chargeof(bunch.species)
    mc2 = massof(bunch.species)
    E_ref = mc2/tilde_m/beta_0
    edge_params = (a, tilde_m, kE, w, w_inv)
    radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing
    params = (radiation_params, kE, beta_0, tilde_m, w, w_inv, a)
    photon_params = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations ? (BeamTracking.exact_elsep!, get_backend(bunch.coords.v), q, mc2, E_ref, kE, w, w_inv) : nothing
    # The potential is taken as zero at the element ends. At the ends of a slice inside the element 
    # (n_slices > 1) it is in general not, so there convert with it as the callbacks do
    if ele.universalparams.first_slice
      kc = push(kc, make_kernel_call(BeamTracking.bmad_to_mad!, (beta_0, tilde_m, 0)))
    else
      kc = push(kc, make_kernel_call(BeamTracking.callback_electric!, (0, 0, beta_0, tilde_m, kE, w, Val{true}())))
    end
    kc = push(kc, integration_launcher(BeamTracking.elsep_rad!, params, photon_params, edge_params, tm, ele.universalparams))
    if ele.universalparams.last_slice
      kc = push(kc, make_kernel_call(BeamTracking.mad_to_bmad!, (beta_0, tilde_m, 0)))
    else
      kc = push(kc, make_kernel_call(BeamTracking.callback_electric!, (0, 0, beta_0, tilde_m, kE, w, Val{false}())))
    end
    kc = push_transforms_out(kc, make_kernel_call(BeamTracking.callback_electric!, (beta_0, tilde_m, kE, w, Val{false}())))
    kc = push_transforms_in(kc, make_kernel_call(BeamTracking.callback_electric!, (beta_0, tilde_m, kE, w, Val{true}())))
    return kc
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
    knl, ksl = get_integrated_strengths(ele.bmultipoleparams, 0, p_over_q_ref, ele.ref.field_scale)
    if isnothing(bunch.coords.q)
      params = (mm, knl, ksl, 0, 0, 0)
    else
      params = (mm, knl, ksl, a, 0, tilde_m)
    end
    return push(kc, make_kernel_call(BeamTracking.integrate_thin!, params))
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
    kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)
    ks[1] == 0 || error("Exact is not currently compatible in an element with both BendParams and a skew dipole multipole")
    Kn0 = kn[1]
    edge_params = (a, tilde_m, beta_0, Kn0, w, w_inv, ele.bendparams.e1, ele.bendparams.e2)
    params = (radiation_params, tilde_m, beta_0, a, g, w, w_inv, Kn0, mm, kn, ks)
    photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, g, w, w_inv, mm, kn, ks) : nothing
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, photon_params, edge_params, tm, ele.universalparams))
  end

  mm = ele.bmultipoleparams.order
  length(mm) == 1 || error("Exact only supports a pure solenoid or pure dipole in a nonzero length element. Try DriftKick, SolenoidKick, MatrixKick or BendKick instead")
  kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)
  photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, 0, nothing, nothing, mm, kn, ks) : nothing

  if first(mm) == 0
    # Pure solenoid
    edge_params = (a, tilde_m, kn[1], nothing, nothing, nothing, nothing, nothing, nothing)
    params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, kn[1], mm, kn, ks)
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