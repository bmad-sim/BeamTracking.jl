
# preparing: the statistics of the Bunch are a PreParameter, computed at the call
@inline function ibs_kick(kc, bunch, ibsparams, bendparams, L, preparing)
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

  ibs_damping_on = (ibsparams.ibs_model != IBSModel.GaussianFluctuationsOnly)
  ibs_fluctuations_on = (ibsparams.ibs_model != IBSModel.GaussianDampingOnly)
  if preparing
    stats = PreParameter((b, L) -> ibs_stats(b, L, Val{true}()), (L,))
  else
    stats = ibs_stats(bunch, L, Val{false}())
    isnothing(stats) && return kc
  end
  params = (get_backend(bunch.coords.v), Val{ibs_damping_on}(), Val{ibs_fluctuations_on}(), stats, g, w, w_inv, L)
  return push(kc, make_kernel_call(BeamTracking.ibs_damping_and_diffusion!, params))
end

# The statistics of the live particles of the Bunch for the IBS kick (nothing without live particles,
# unless prepared)
function ibs_stats(bunch, L, ::Val{prepared}) where {prepared}
  log_c_light = log(C_LIGHT)
  log_m = log(massof(bunch.species)) - 2*log_c_light
  log_q = log(abs(chargeof(bunch.species)))
  log_k = log_m + 2*log_q - log(4*pi*EPS_0)
  # Lost particles do not count (their coordinates are frozen where they were lost)
  T = eltype(bunch.coords.v)
  weights = (bunch.coords.state .== BeamTracking.STATE_ALIVE) .*
            (isnothing(bunch.coords.weight) ? one(T) : bunch.coords.weight isa Number ? T(bunch.coords.weight) : bunch.coords.weight)
  W = sum(weights)
  # No live particle: nothing to kick (when prepared, any finite statistics: the kernel only changes
  # live particles)
  if iszero(W)
    prepared || return nothing
    weights = one.(weights)
    W = sum(weights)
  end
  log_N = log(W)
  tilde_m, gamsqr_0, _ = drift_params(bunch.species, bunch.p_over_q_ref)
  log_p0 = log_m + log_c_light - log(tilde_m)
  gamma_0 = sqrt(gamsqr_0)
  backend = get_backend(bunch.coords.v)

  means, sigma = mean_and_cov(bunch.coords.v, weights, backend)
  sigma = Symmetric(sigma)
  sigma_inv = inv(sigma)
  log_b_min = log(4) + log_k - log(sigma[2,2] + sigma[4,4] + sigma[6,6]/gamsqr_0) - 2*log_p0
  log_b_max = log(minimum((sigma[1,1], sigma[3,3], sigma[5,5]*gamsqr_0)))/2
  L_C = log_b_max - log_b_min

  C = SA[sigma_inv[2,2]         sigma_inv[2,4]         sigma_inv[2,6]*gamma_0;
         sigma_inv[4,2]         sigma_inv[4,4]         sigma_inv[4,6]*gamma_0;
         sigma_inv[6,2]*gamma_0 sigma_inv[6,4]*gamma_0 sigma_inv[6,6]*gamsqr_0]
  lambdas, vectors = eigen(Symmetric(C))
  P = vectors'
  integrals = ibs_integrals(lambdas...)

  log_minus_b_coeff = log_N + 2*log_k + log(L_C/pi^2) - log_m - 3*log_p0 - logdet(sigma)/2
  b_coeff = -exp(log_minus_b_coeff)
  d_coeff = -b_coeff/2

  d_xx = d_coeff/gamma_0*((1-P[1,1]^2)*integrals[1] + (1-P[2,1]^2)*integrals[2] + (1-P[3,1]^2)*integrals[3])
  d_xy = d_coeff/gamma_0*(-P[1,1]*P[1,2]*integrals[1] - P[2,1]*P[2,2]*integrals[2] - P[3,1]*P[3,2]*integrals[3])
  d_xz = d_coeff*(-P[1,1]*P[1,3]*integrals[1] - P[2,1]*P[2,3]*integrals[2] - P[3,1]*P[3,3]*integrals[3])
  d_yy = d_coeff/gamma_0*((1-P[1,2]^2)*integrals[1] + (1-P[2,2]^2)*integrals[2] + (1-P[3,2]^2)*integrals[3])
  d_yz = d_coeff*(-P[1,2]*P[1,3]*integrals[1] - P[2,2]*P[2,3]*integrals[2] - P[3,2]*P[3,3]*integrals[3])
  d_zz = d_coeff*gamma_0*((1-P[1,3]^2)*integrals[1] + (1-P[2,3]^2)*integrals[2] + (1-P[3,3]^2)*integrals[3])

  diffusion_mat = SA[d_xx d_xy d_xz;
                     d_xy d_yy d_yz;
                     d_xz d_yz d_zz]
  diffusion_lambdas, diffusion_vectors = eigen(Symmetric(diffusion_mat))
  diffusion_P = diffusion_vectors'

  sigma_inv_t = (sigma_inv[1,1], sigma_inv[1,2], sigma_inv[1,3], sigma_inv[1,4], sigma_inv[1,5], sigma_inv[1,6],
                                 sigma_inv[2,2], sigma_inv[2,3], sigma_inv[2,4], sigma_inv[2,5], sigma_inv[2,6],
                                                 sigma_inv[3,3], sigma_inv[3,4], sigma_inv[3,5], sigma_inv[3,6],
                                                                 sigma_inv[4,4], sigma_inv[4,5], sigma_inv[4,6],
                                                                                 sigma_inv[5,5], sigma_inv[5,6],
                                                                                                 sigma_inv[6,6])
  return (tilde_m, gamma_0, b_coeff, integrals, diffusion_lambdas, diffusion_P, P, sigma_inv_t, means)
end

#---------------------------------------------------------------------------------------------------

# s: position of the face (entrance at s = 0, exit at s = L, or a slice end) within the element
@inline function alignment(kc, alignmentparams, bendparams, L, s, entering)
  if isnothing(alignmentparams); return kc; end

  x_off = alignmentparams.x_offset
  y_off = alignmentparams.y_offset
  z_off = alignmentparams.z_offset
  x_rot = alignmentparams.x_rot
  y_rot = alignmentparams.y_rot
  tilt  = alignmentparams.tilt

  ele_orient = 1   ## Future work: Need to extend this for reversed elements.

  if !isnothing(bendparams) # && (bendparams.g_ref != 0 || bendparams.tilt_ref != 0)
    if entering
      mid_r, mid_q, st, ct = BeamTracking.coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, L)
      kc = push(kc, make_kernel_call(BeamTracking.track_coord_bend_transform_at_s!, (mid_r, mid_q, st, ct, bendparams.g_ref, L, s, Val{true}())))
      kc = push_transforms_in(kc, make_kernel_call(
          BeamTracking.callback_coord_bend_transform_at_s!, (mid_r, mid_q, st, ct,  bendparams.g_ref, L, Val{true}())
        )
      )
      return kc
    else
      mid_r, mid_q, st, ct = BeamTracking.coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, L)
      kc = push(kc, make_kernel_call(BeamTracking.track_coord_bend_transform_at_s!, (mid_r, mid_q, st, ct, bendparams.g_ref, L, s, Val{false}())))
      kc = push_transforms_out(kc, make_kernel_call(
          BeamTracking.callback_coord_bend_transform_at_s!, (mid_r, mid_q, st, ct,  bendparams.g_ref, L, Val{false}())
        )
      )
      return kc
    end
  else
    if entering
      kc = push(kc, make_kernel_call(BeamTracking.track_alignment_straight_at_s!, (x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, s, Val{true}())))
      kc = push_transforms_in(kc, make_kernel_call(
          BeamTracking.callback_alignment_straight_at_s!, (x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, Val{true}())
        )
      )
      return kc
    else
      kc = push(kc, make_kernel_call(BeamTracking.track_alignment_straight_at_s!, (x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, s, Val{false}())))
      kc = push_transforms_out(kc, make_kernel_call(
          BeamTracking.callback_alignment_straight_at_s!, (x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L, Val{false}())
        )
      )
      return kc
    end
  end
end

#---------------------------------------------------------------------------------------------------

@inline function aperture(kc, apertureparams, entering)
  x1 = apertureparams.x1_limit
  x2 = apertureparams.x2_limit
  y1 = apertureparams.y1_limit
  y2 = apertureparams.y2_limit

  if entering && apertureparams.aperture_at == ApertureAt.Exit
    return kc
  elseif !entering && apertureparams.aperture_at == ApertureAt.Entrance
    return kc
  elseif apertureparams.aperture_shape == ApertureShape.Elliptical
    if any(isinf, (x1, x2, y1, y2))
      error("Invalid ApertureParams limits for elliptical aperture: check if all limits have been set")
    end
    return push(kc, make_kernel_call(BeamTracking.track_aperture_elliptical!, (x1, x2, y1, y2)))
  else  
    return push(kc, make_kernel_call(BeamTracking.track_aperture_rectangular!, (x1, x2, y1, y2)))
  end
end
#---------------------------------------------------------------------------------------------------
