#---------------------------------------------------------------------------------------------------

@inline function alignment(tm, bunch, alignmentparams, bendparams, L, entering)
  if !isactive(alignmentparams); return nothing; end

  x_off = alignmentparams.x_offset
  y_off = alignmentparams.y_offset
  z_off = alignmentparams.z_offset
  x_rot = alignmentparams.x_rot
  y_rot = alignmentparams.y_rot
  tilt  = alignmentparams.tilt

  ele_orient = 1   ## Future work: Need to extend this for reversed elements.

  #

  if isactive(bendparams) && (bendparams.g_ref != 0 || bendparams.tilt_ref != 0)
    if entering
      dr, q = BeamTracking.coord_alignment_bend_entering(x_off, y_off, z_off, 
                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, L)
      return KernelCall(BeamTracking.track_coord_transform!, (dr, q))
    else
      dr, q = BeamTracking.coord_alignment_bend_exiting(x_off, y_off, z_off, 
                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, L)
      return KernelCall(BeamTracking.track_coord_transform!, (dr, q))
    end

  else
    if entering
      return KernelCall(BeamTracking.track_alignment_straight_entering!, (x_off, y_off, z_off, 
                                                     x_rot, y_rot, tilt, ele_orient, L))
    else
      return KernelCall(BeamTracking.track_alignment_straight_exiting!, (x_off, y_off, z_off, 
                                                     x_rot, y_rot, tilt, ele_orient, L))
    end
  end
end

#---------------------------------------------------------------------------------------------------

@inline function aperture(tm, bunch, apertureparams, entering)
  x1 = apertureparams.x1_limit
  x2 = apertureparams.x2_limit
  y1 = apertureparams.y1_limit
  y2 = apertureparams.y2_limit

  if entering && apertureparams.aperture_at == ApertureAt.Exit
      return KernelCall()
  elseif !entering && apertureparams.aperture_at == ApertureAt.Entrance
      return KernelCall()
  elseif apertureparams.aperture_shape == ApertureShape.Elliptical
    if any(isinf, (x1, x2, y1, y2))
      error("Invalid ApertureParams limits for elliptical aperture: check if all limits have been set")
    end
    return KernelCall(BeamTracking.track_aperture_elliptical!, (x1, x2, y1, y2))
  else  
    return KernelCall(BeamTracking.track_aperture_rectangular!, (x1, x2, y1, y2))
  end
end

#---------------------------------------------------------------------------------------------------

@inline function RFcavity(tm, bunch, bmultipoleparams, rfparams, beamlineparams, L)
  if !isactive(bmultipoleparams)
    return pure_rf(tm, bunch, rfparams, beamlineparams, L)
  else
    return bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, beamlineparams, L)
  end
end

#---------------------------------------------------------------------------------------------------

@inline function ibs_kick(tm, bunch, bendparams, L)
  p_over_q_ref = bunch.p_over_q_ref
  if !isnothing(bendparams)
    g = bendparams.g_ref
    tilt = bendparams.tilt_ref
    w = rot_quaternion(0, 0, -tilt)
    w_inv = inv_rot_quaternion(0, 0, -tilt)
  else
    g = 0
    w = nothing
    w_inv = nothing
  end

  log_e_charge = log(E_CHARGE)
  log_c_light = log(C_LIGHT)
  log_m = log(massof(bunch.species)) - 2*log_c_light + log_e_charge
  log_q = log(abs(chargeof(bunch.species))) + log_e_charge
  log_k = log_m + 2*log_q - log(4*pi*EPS_0)
  log_N = log(tm.ibs_num_particles)
  tilde_m, gamsqr_0, _ = BeamTracking.drift_params(bunch.species, p_over_q_ref)
  log_p0 = log_m + log_c_light - log(tilde_m)
  gamma_0 = sqrt(gamsqr_0)

  means, sigma = mean_and_cov(bunch.coords)
  sigma = Symmetric(sigma)
  log_b_min = log(4) + log_k - log(sigma[2,2] + sigma[4,4] + sigma[6,6]/gamsqr_0) - 2*log_p0
  log_b_max = log(minimum((sigma[1,1], sigma[3,3], sigma[5,5]*gamsqr_0)))/2
  L_C = log_b_max - log_b_min
  sigma_inv = inv(sigma)
  C = SA[sigma_inv[2,2]         sigma_inv[2,4]         sigma_inv[2,6]*gamma_0;
         sigma_inv[4,2]         sigma_inv[4,4]         sigma_inv[4,6]*gamma_0;
         sigma_inv[6,2]*gamma_0 sigma_inv[6,4]*gamma_0 sigma_inv[6,6]*gamsqr_0]
  lambdas, vectors = eigen(Symmetric(C))
  P = vectors'
  integrals = ibs_integrals(lambdas...)
  #println(lambdas)
  #println(integrals)
  
  log_minus_b_coeff = log_N + 2*log_k + log(L_C/pi^2) - log_m - 3*log_p0 - logdet(sigma)/2
  b_coeff = -exp(log_minus_b_coeff)
  d_coeff = -b_coeff/2
  #println(d_coeff)

  I_XX = (1-P[1,1]^2)*integrals[1] + (1-P[2,1]^2)*integrals[2] + (1-P[3,1]^2)*integrals[3]
  I_YY = (1-P[1,2]^2)*integrals[1] + (1-P[2,2]^2)*integrals[2] + (1-P[3,2]^2)*integrals[3]
  I_ZZ = (1-P[1,3]^2)*integrals[1] + (1-P[2,3]^2)*integrals[2] + (1-P[3,3]^2)*integrals[3]

  D_XX = d_coeff*I_XX
  D_YY = d_coeff*I_YY
  D_ZZ = d_coeff*I_ZZ

  d_xx = D_XX/gamma_0
  d_yy = D_YY/gamma_0
  d_zz = D_ZZ*gamma_0
  diffusion = (d_xx, d_yy, d_zz)
  #println(diffusion)
  
  params = (tilde_m, gamma_0, b_coeff, integrals, diffusion, P, sigma_inv, means, g, w, w_inv, L)
  return KernelCall(BeamTracking.ibs_damping_and_diffusion!, params)
end