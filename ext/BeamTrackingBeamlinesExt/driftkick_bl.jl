@inline function fringe_body_fringe(tm::DriftKick, kc, p_over_q_ref, bunch, ele)
  L = ele.universalparams.L
  isnothing(ele.emultipoleparams) || error("DriftKick is not compatible with elements containing EMultipoleParams. For an electric dipole, try Exact instead")
  isnothing(ele.bendparams) || error("DriftKick is not compatible with elements containing BendParams. Try BendKick instead")
  isnothing(ele.mapparams) || error("DriftKick is not compatible with elements containing MapParams. Try Exact instead")
  isnothing(ele.patchparams) || error("DriftKick is not compatible with elements containing PatchParams. Try Exact instead")
  isnothing(ele.fourpotentialparams) || error("DriftKick is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.emfieldparams) || error("DriftKick is not compatible with elements containing EMFieldParams. Try RK4 instead")
  L != 0 || error("DriftKick is not compatible with zero length elements. Try Exact instead")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)

  if isnothing(ele.bmultipoleparams) && isnothing(ele.rfparams)
    # Pure drift, use order 2 (exact solution)
    tm_new = @set tm.order = 2
    params = (beta_0, gamsqr_0, tilde_m)
    return push(kc, integration_launcher(BeamTracking.exact_drift!, params, nothing, nothing, tm_new, ele.universalparams))
  end

  q = chargeof(bunch.species)
  mc2 = massof(bunch.species)
  a = gyromagnetic_anomaly(bunch.species)
  E_ref = mc2/tilde_m/beta_0
  radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing

  if !isnothing(ele.rfparams)
    omega = rf_omega_calc(ele.rfparams, ele.beamlineparams)
    t_ref = (rf_phi0_calc(ele.rfparams, ele.beamlineparams.beamline.species_ref) - pi/2)/omega
    E0_normalized = ele.rfparams.voltage/L/p_over_q_ref
  end

  if isnothing(ele.bmultipoleparams)
    # Cavity with no BMultipoles
    params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, omega, t_ref, E0_normalized, nothing, Val{false}(), SA[], SA[], SA[])
    return push(kc, integration_launcher(BeamTracking.cavity!, params, nothing, nothing, tm, ele.universalparams))
  end

  mm = ele.bmultipoleparams.order
  first(mm) != 0 || error("DriftKick is not compatible in element with solenoidal multipole component. Try SolenoidKick instead")
  kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref)
  Kn0::Union{Nothing,eltype(kn)} = nothing
  Ks0::Union{Nothing,eltype(ks)} = nothing
  Kn1::Union{Nothing,eltype(kn)} = nothing
  Ks1::Union{Nothing,eltype(ks)} = nothing
  Kn2::Union{Nothing,eltype(kn)} = nothing
  Ks2::Union{Nothing,eltype(ks)} = nothing
  for j in 1:length(mm)
    if mm[j] == 1
      Kn0 = kn[j]
      Ks0 = ks[j]
    elseif mm[j] == 2
      Kn1 = kn[j]
      Ks1 = ks[j]
    elseif mm[j] == 3
      Kn2 = kn[j]
      Ks2 = ks[j]
    end
  end
  edge_params = (a, tilde_m, nothing, Kn0, Ks0, Kn1, Ks1, Kn2, Ks2)
  photon_params = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, 0, nothing, nothing, mm, kn, ks) : nothing

  if isnothing(ele.rfparams)
    # Just BMultipoles (not a cavity)
    params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks)
    return push(kc, integration_launcher(BeamTracking.dkd_multipole!, params, photon_params, edge_params, tm, ele.universalparams))
  else
    # Cavity with BMultipoles
    params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, omega, t_ref, E0_normalized, nothing, Val{false}(), mm, kn, ks)
    return push(kc, integration_launcher(BeamTracking.cavity!, params, photon_params, edge_params, tm, ele.universalparams))
  end
end