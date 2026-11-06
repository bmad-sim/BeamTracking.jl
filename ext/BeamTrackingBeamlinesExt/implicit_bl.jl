@inline function fringe_body_fringe(tm::Implicit, kc, p_over_q_ref, bunch, ele)
  L = ele.universalparams.L
  isnothing(ele.bmultipoleparams) || error("Implicit is not compatible with elements containing BMultipoleParams. Include the multipole fields in the FourPotentialParams instead")
  isnothing(ele.emultipoleparams) || error("Implicit is not compatible with elements containing EMultipoleParams. Include the multipole fields in the FourPotentialParams instead")
  isnothing(ele.rfparams) || error("Implicit is not compatible with elements containing RFParams. Try SolenoidKick or DriftKick instead")
  isnothing(ele.patchparams) || error("Implicit is not compatible with elements containing PatchParams. Try Exact instead")
  isnothing(ele.mapparams) || error("Implicit is not compatible with elements containing MapParams. Try Exact instead")
  isnothing(ele.emfieldparams) || error("Implicit is not compatible with elements containing EMFieldParams. Try RK4 instead")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)

  if isnothing(ele.fourpotentialparams)
    # Pure drift, use order 2 (exact solution)
    isnothing(ele.bendparams) || error("Implicit is not compatible with elements containing BendParams but no FourPotentialParams. Try BendKick instead")
    tm_new = @set tm.order = 2
    params = (beta_0, gamsqr_0, tilde_m)
    return push(kc, integration_launcher(BeamTracking.exact_drift!, params, nothing, nothing, tm_new, ele.universalparams))
  end

  if isnothing(ele.bendparams)
    g = zero(L)
    w = nothing
    w_inv = nothing
  else
    bp = ele.bendparams
    bp.e1 == 0 && bp.e2 == 0 || error("Implicit is not compatible with elements with nonzero bend edge angles e1 or e2")
    g = bp.g_ref
    ntilt = -bp.tilt_ref
    w = ntilt == 0 ? nothing : rot_quaternion(0, 0, ntilt)
    w_inv = ntilt == 0 ? nothing : inv_rot_quaternion(0, 0, ntilt)
  end

  fpp = ele.fourpotentialparams
  # Normalized potentials are referenced to the Beamline p_over_q_ref (see field_scale)
  !fpp.four_potential_normalized || ele.ref.field_scale == 1 || error("Implicit is not yet compatible with normalized FourPotentialParams when the Bunch species has the opposite charge sign of the reference species")
  potential_and_jac = fpp.four_potential
  potential_params = fpp.four_potential_params
  normalized = fpp.four_potential_normalized ? Val{true}() : Val{false}()
  use_newton = tm.use_newton ? Val{true}() : Val{false}()

  q = chargeof(bunch.species)
  mc2 = massof(bunch.species)
  a = gyromagnetic_anomaly(bunch.species)
  E_ref = mc2/tilde_m/beta_0
  radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing
  # The potential is evaluated at s from the element entrance (s0 is nonzero for a slice of the element)
  s0 = ele.universalparams.s0
  params = (radiation_params, beta_0, tilde_m, a, g, w, w_inv, potential_and_jac, potential_params, p_over_q_ref, normalized, use_newton, s0)
  photon_params = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations ? (BeamTracking.implicit_integrator!, get_backend(bunch.coords.v), q, mc2, E_ref, g, w, w_inv, potential_and_jac, potential_params, p_over_q_ref, normalized, s0) : nothing

  # The scalar potential is taken as zero at the element ends. At the ends of a slice inside the
  # element (n_slices > 1) it is in general not, so there convert with it as the callbacks do
  if ele.universalparams.first_slice
    kc = push(kc, make_kernel_call(BeamTracking.bmad_to_mad!, (beta_0, tilde_m, 0)))
  else
    kc = push(kc, make_kernel_call(BeamTracking.callback_implicit!, (s0, 0, beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, Val{true}())))
  end
  kc = push(kc, integration_launcher(BeamTracking.implicit_integrator!, params, photon_params, nothing, tm, ele.universalparams))
  if ele.universalparams.last_slice
    kc = push(kc, make_kernel_call(BeamTracking.mad_to_bmad!, (beta_0, tilde_m, 0)))
  else
    kc = push(kc, make_kernel_call(BeamTracking.callback_implicit!, (s0 + ele.universalparams.L_slice, 0, beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, Val{false}())))
  end
  kc = push_transforms_out(kc, make_kernel_call(BeamTracking.callback_implicit!, (beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, Val{false}())))
  kc = push_transforms_in(kc, make_kernel_call(BeamTracking.callback_implicit!, (beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, Val{true}())))
  return kc
end