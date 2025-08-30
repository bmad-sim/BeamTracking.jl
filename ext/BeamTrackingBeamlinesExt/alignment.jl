@inline function alignment(tm, bunch, alignmentparams, bendparams, L, entering::Bool)
#  if !isactive(alignmentparams); return nothing; end
#
#  x_off = alignmentparams.x_off
#  y_off = alignmentparams.y_off
#  z_off = alignmentparams.z_off
#  x_rot = alignmentparams.x_rot
#  y_rot = alignmentparams.y_rot
#  tilt  = alignmentparams.tilt
#  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
#  ele_orient = 1   # Need to extend this to reversed elements.
#
#  if isactive(bendparams) && bendparams.g_ref != 0
#    return KernelCall(track_alignment_bend!, (entering, x_off, y_off, z_off, 
#                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, tilde_m, beta_0, L))
#  else
#    return KernelCall(track_alignment_straight!, (entering, x_off, y_off, z_off, 
#                                         x_rot, y_rot, tilt, ele_orient, tilde_m, beta_0, L))
#  end
end