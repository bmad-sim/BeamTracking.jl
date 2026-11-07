@inline function fringe_body_fringe(tm::MatrixKick, kc, bunch, ele)
  L = ele.universalparams.L
  p_over_q_ref = ele.ref.p_over_q_ref
  isnothing(ele.emultipoleparams) || error("MatrixKick is not compatible with elements containing EMultipoleParams. For an electric dipole, try Exact instead")
  isnothing(ele.bendparams) || error("MatrixKick is not compatible with elements containing BendParams. Try BendKick instead")
  isnothing(ele.mapparams) || error("MatrixKick is not compatible with elements containing MapParams. Try Exact instead")
  isnothing(ele.patchparams) || error("MatrixKick is not compatible with elements containing PatchParams. Try Exact instead")
  isnothing(ele.fourpotentialparams) || error("MatrixKick is not compatible with elements containing FourPotentialParams. Try Implicit instead")
  isnothing(ele.rfparams) || error("MatrixKick is not compatible with elements containing RFParams. Try SolenoidKick or DriftKick instead")
  isnothing(ele.emfieldparams) || error("MatrixKick is not compatible with elements containing EMFieldParams. Try RK4 instead")
  L != 0 || error("MatrixKick is not compatible with zero length elements. Try Exact instead")

  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)

  if isnothing(ele.bmultipoleparams)
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

  mm = ele.bmultipoleparams.order
  first(mm) != 0 || error("MatrixKick is not compatible in element with solenoidal multipole component. Try SolenoidKick instead")
  kn, ks = get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)
  Kn0::Union{Nothing,eltype(kn)} = nothing
  Ks0::Union{Nothing,eltype(ks)} = nothing
  # The type of a TimeDependentParam changes with every operation (zero, sqrt, ...)
  TN = eltype(kn) <: TimeDependentParam ? TimeDependentParam : eltype(kn)
  TS = eltype(ks) <: TimeDependentParam ? TimeDependentParam : eltype(ks)
  Kn1::TN = zero(kn[1])
  Ks1::TS = zero(ks[1])
  k1::TN = zero(kn[1])
  tilt1::TN = zero(kn[1])
  Kn2::Union{Nothing,eltype(kn)} = nothing
  Ks2::Union{Nothing,eltype(ks)} = nothing
  for j in 1:length(mm)
    if mm[j] == 1
      Kn0 = kn[j]
      Ks0 = ks[j]
    elseif mm[j] == 2
      Kn1 = kn[j]
      Ks1 = ks[j]
      k1 = sqrt(Kn1^2 + Ks1^2)
      tilt1 = atan2(Ks1, Kn1) / 2
    elseif mm[j] == 3
      Kn2 = kn[j]
      Ks2 = ks[j]
    end
  end
  if tilt1 == 0
    w1 = nothing
    w1_inv = nothing
  else
    w1 = rot_quaternion(0, 0, tilt1)
    w1_inv = inv_rot_quaternion(0, 0, tilt1)
  end
  edge_params = (a, tilde_m, nothing, Kn0, Ks0, Kn1, Ks1, Kn2, Ks2)
  params = (radiation_params, beta_0, gamsqr_0, tilde_m, a, w1, w1_inv, k1, mm, kn, ks)
  photon_params = isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations ? (get_backend(bunch.coords.v), q, mc2, E_ref, 0, nothing, nothing, mm, kn, ks) : nothing
  
  return push(kc, integration_launcher(BeamTracking.mkm_quadrupole!, params, photon_params, edge_params, tm, ele.universalparams))
end