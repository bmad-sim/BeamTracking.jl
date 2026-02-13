@inline function alignment(tm::Exact, bunch, alignmentparams, bendparams, L, entering)
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


@inline function aperture(tm::Exact, bunch, apertureparams, entering)
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


@inline function pure_patch(tm::Exact, bunch, patchparams, L) 
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  winv = inv_rot_quaternion(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
  return KernelCall(BeamTracking.patch!, (beta_0, gamsqr_0, tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L))
end

@inline function thick_pure_bsolenoid(tm::Exact, bunch, bm0, L)
  Ksol, _ = get_strengths(bm0, L, bunch.p_over_q_ref)
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_solenoid!, (Ksol, beta_0, gamsqr_0, tilde_m, L))
end

@inline function drift(tm::Exact, bunch, L)
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_drift!, (beta_0, gamsqr_0, tilde_m, L))
end

@inline function thick_bend_pure_bdipole(tm::Exact, bunch, bendparams, bm1, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  if tm.fringe_at == Fringe.BothEnds || tm.fringe_at == Fringe.EntranceEnd
    e1 = bendparams.e1
  else
    e1 = 0
  end
  if tm.fringe_at == Fringe.BothEnds || tm.fringe_at == Fringe.ExitEnd
    e2 = bendparams.e2
  else
    e2 = 0
  end
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  Kn0, Ks0 = get_strengths(bm1, L, bunch.p_over_q_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_bend_with_rotation!, (e1, e2, theta, 0, g, Kn0, w, w_inv, tilde_m, beta_0, L))
end

@inline function thick_pure_bdipole(tm::Exact, bunch, bm1, L)
  Kn0, Ks0 = get_strengths(bm1, L, bunch.p_over_q_ref)
  Kn = sqrt(Kn0^2 + Ks0^2)
  tilt = atan2(Ks0, Kn0)
  w = rot_quaternion(0,0,tilt)
  w_inv = inv_rot_quaternion(0,0,tilt)
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_bend_with_rotation!, (0, 0, 0, 0, 0, Kn, w, w_inv, tilde_m, beta_0, L))
end

@inline function thick_bend_no_field(tm::Exact, bunch, bendparams, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  tilde_m, _, beta_0 = BeamTracking.drift_params(bunch.species, bunch.p_over_q_ref)
  return KernelCall(BeamTracking.exact_curved_drift!, (e1, e2, theta, g, w, w_inv, gyromagnetic_anomaly(bunch.species), tilde_m, beta_0, L))
end

@inline pure_map(tm::Exact, bunch, mapparams, L) = KernelCall(BeamTracking.map!, (mapparams.transport_map, mapparams.transport_map_params, L))

