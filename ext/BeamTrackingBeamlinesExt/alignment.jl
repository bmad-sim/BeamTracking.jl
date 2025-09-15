@inline function alignment(tm, bunch, alignmentparams, bendparams, L, entering::Bool)
  if !isactive(alignmentparams); return nothing; end

  x_off = alignmentparams.x_off
  y_off = alignmentparams.y_off
  z_off = alignmentparams.z_off
  x_rot = alignmentparams.x_rot
  y_rot = alignmentparams.y_rot
  tilt  = alignmentparams.tilt

  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  ele_orient = 1   ## Future work: Need to extend this for reversed elements.

  if isactive(bendparams) && bendparams.g_ref != 0
    return KernelCall(track_alignment_bend!, (beta_0, gamsqr_0, tilde_m, entering, x_off, y_off, z_off, 
                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, L))
  else
    return KernelCall(track_alignment_straight!, (beta_0, gamsqr_0, tilde_m, entering, x_off, y_off, z_off, 
                                         x_rot, y_rot, tilt, ele_orient, L))
  end
end