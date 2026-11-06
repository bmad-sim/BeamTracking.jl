@inline function fringe_body_fringe(tm::BendKick, kc, p_over_q_ref, bunch, ele)
  L = ele.universalparams.L
  isnothing(ele.emultipoleparams) || error("BendKick is not compatible with elements containing EMultipoleParams. For an electric dipole, try Exact instead")
  isnothing(ele.mapparams) || error("BendKick is not compatible with elements containing MapParams. Try Exact instead")
  isnothing(ele.patchparams) || error("BendKick is not compatible with elements containing PatchParams. Try Exact instead")
  isnothing(ele.fourpotentialparams) || error("BendKick is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.rfparams) || error("BendKick is not compatible with elements containing RFParams")
  isnothing(ele.emfieldparams) || error("BendKick is not compatible with elements containing EMFieldParams. Try RK4 instead")
  L != 0 || error("BendKick is not compatible with zero length elements. Try Exact instead")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)

  if isnothing(ele.bendparams) && isnothing(ele.bmultipoleparams)
    # Pure drift, use order 2 (exact solution)
    tm_new = @set tm.order = 2
    params = (beta_0, gamsqr_0, tilde_m)
    return push(kc, integration_launcher(BeamTracking.exact_drift!, params, nothing, nothing, tm_new, ele.universalparams))
  end

  a = gyromagnetic_anomaly(bunch.species)

  if isnothing(ele.bmultipoleparams)
    # Pure bend, use order 2 (exact solution)
    tm_new = @set tm.order = 2
    g = ele.bendparams.g_ref
    ntilt = -ele.bendparams.tilt_ref
    w = ntilt == 0 ? nothing : rot_quaternion(0, 0, ntilt)
    w_inv = ntilt == 0 ? nothing : inv_rot_quaternion(0, 0, ntilt)
    params = (nothing, tilde_m, beta_0, a, g, w, w_inv, 0, SA[0], SA[0], SA[0])
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, nothing, nothing, tm_new, ele.universalparams))
  end

  mm = ele.bmultipoleparams.order
  kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)
  first(mm) != 0 || error("BendKick is not compatible in element with solenoidal multipole component. Try SolenoidKick instead")
  has_dipole = first(mm) == 1

  q = chargeof(bunch.species)
  mc2 = massof(bunch.species)
  E_ref = mc2/tilde_m/beta_0
  radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing
  fluctuations = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations
  # Pure bend, use order 2 (exact solution). Not with radiation, which is applied in the
  # integrator steps, so the order sets how accurately it is integrated
  tm_new = (has_dipole && length(mm) == 1 && isnothing(radiation_params) && !fluctuations) ? (@set tm.order = 2) : tm

  if isnothing(ele.bendparams)
    # Use straight magnet fringe field model AND obtain "tilt_ref" from tilt0, Kn0, Ks0:
    Kn0 = has_dipole ? sqrt(kn[1]^2 + ks[1]^2) : zero(kn[1]) # In body frame dipole strength
    kn0 = has_dipole ? kn[1] : zero(kn[1]) # normal/skew strengths for straight magnet fringe calc
    ks0 = has_dipole ? ks[1] : zero(ks[1])
    tilt0 = has_dipole ? atan(ks[1], kn[1]) : zero(kn[1])
    Kn1::Union{Nothing,eltype(kn)} = nothing
    Ks1::Union{Nothing,eltype(ks)} = nothing
    Kn2::Union{Nothing,eltype(kn)} = nothing
    Ks2::Union{Nothing,eltype(ks)} = nothing
    for j in 1:length(mm)
      if mm[j] == 2
        Kn1 = kn[j]
        Ks1 = ks[j]
      elseif mm[j] == 3
        Kn2 = kn[j]
        Ks2 = ks[j]
      end
    end
    w = tilt0 == 0 ? nothing : rot_quaternion(0, 0, tilt0)
    w_inv = tilt0 == 0 ? nothing : inv_rot_quaternion(0, 0, tilt0)
    edge_params = (a, tilde_m, nothing, kn0, ks0, Kn1, Ks1, Kn2, Ks2)
    # Strengths in the dipole frame (rotated by w), rotated as get_strengths does for a tilt
    knr = @. kn*cos(mm*tilt0) + ks*sin(mm*tilt0)
    ksr = @. -kn*sin(mm*tilt0) + ks*cos(mm*tilt0)
    params = (radiation_params, tilde_m, beta_0, a, 0, w, w_inv, Kn0, mm, knr, ksr)
    photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, 0, nothing, nothing, mm, kn, ks) : nothing
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, photon_params, edge_params, tm_new, ele.universalparams))
  else
    g = ele.bendparams.g_ref
    ntilt = -ele.bendparams.tilt_ref
    !has_dipole || ks[1] == 0 || error("BendKick is not currently compatible in an element with both BendParams and a skew dipole multipole")
    has_dipole && length(mm) == 1 || @warn "Straight multipoles are being used in a curved reference system. Maxwell's equations in free space are not satisfied." maxlog=1
    Kn0 = has_dipole ? kn[1] : zero(kn[1])
    w = ntilt == 0 ? nothing : rot_quaternion(0, 0, ntilt)
    w_inv = ntilt == 0 ? nothing : inv_rot_quaternion(0, 0, ntilt)
    edge_params = (a, tilde_m, beta_0, Kn0, w, w_inv, ele.bendparams.e1, ele.bendparams.e2)
    params = (radiation_params, tilde_m, beta_0, a, g, w, w_inv, Kn0, mm, kn, ks)
    photon_params = fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, g, w, w_inv, mm, kn, ks) : nothing
    return push(kc, integration_launcher(BeamTracking.bkb_multipole!, params, photon_params, edge_params, tm_new, ele.universalparams))
  end
end