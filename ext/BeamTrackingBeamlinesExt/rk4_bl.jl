@inline function fringe_body_fringe(tm::RK4, kc, bunch, ele)
  L = ele.universalparams.L
  p_over_q_ref = ele.ref.p_over_q_ref
  isnothing(bunch.coords.q) || error("RK4 is not compatible with spin tracking")
  L > 0 || error("RK4 is not compatible with zero or negative length elements. Try Exact instead")
  isnothing(ele.patchparams) || error("RK4 is not compatible with elements containing PatchParams. Try Exact instead")
  isnothing(ele.rfparams) || error("RK4 is not compatible with elements containing RFParams. Try SolenoidKick or DriftKick instead")
  isnothing(ele.mapparams) || error("RK4 is not compatible with elements containing MapParams. Try Exact instead")
  isnothing(ele.fourpotentialparams) || error("RK4 is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.emultipoleparams) || error("RK4 is not compatible with elements containing EMultipoleParams. For an electric dipole, try Exact instead")

  if isnothing(ele.bendparams)
    gx = zero(L)
    gy = zero(L)
  else
    bp = ele.bendparams
    bp.edge1_int == 0 && bp.edge2_int == 0 || error("RK4 is not compatible with elements with nonzero edge1_int or edge2_int")
    bp.e1 == 0 && bp.e2 == 0 || error("RK4 is not compatible with elements with nonzero bend edge angles e1 or e2. Try BendKick or Exact instead")
    gx = bp.g_ref * cos(bp.tilt_ref)
    gy = bp.g_ref * sin(bp.tilt_ref)
  end

  if isnothing(ele.bmultipoleparams)
    multipole_functions = ()
    multipole_parameters = ()
    multipole_normalized = ()
  else
    mm = ele.bmultipoleparams.order
    bn, bs = get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)
    multipole_functions = (BeamTracking.multipole_field,)
    multipole_parameters = ((mm, bn, bs),)
    multipole_normalized = (Val{true}(),)
  end

  ep = ele.emfieldparams
  if isnothing(ep) || isnothing(ep.em_field)
    custom_functions = ()
    custom_parameters = ()
    custom_normalized = ()
  else
    custom_functions = (ep.em_field,)
    custom_parameters = (ep.em_field_params,)
    # Normalized fields are normalized by the Beamline p_over_q_ref, i.e. the tracking p_over_q_ref
    # times field_scale: for field_scale = -1, pass the factor -1 (see normalized_field_at)
    if !ep.em_field_normalized
      custom_normalized = (Val{false}(),)
    elseif ele.ref.field_scale == 1
      custom_normalized = (Val{true}(),)
    else
      custom_normalized = (-1.0,)
    end
  end

  species = bunch.species
  tilde_m, _, beta_0 = drift_params(species, p_over_q_ref)
  charge = chargeof(species)
  p0c = R_to_pc(species, p_over_q_ref)
  mc2 = massof(species)
  n_steps = ele.universalparams.n_steps
  ds_step = ele.universalparams.ds_step

  # The fields are evaluated at s from the element entrance (s0 is nonzero for a slice of the element)
  params = (beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps, gx, gy,
            (multipole_functions..., custom_functions...),
            (multipole_parameters..., custom_parameters...),
            (multipole_normalized..., custom_normalized...), ele.universalparams.s0)
  return push(kc, make_kernel_call(BeamTracking.rk4_kernel!, params))
end