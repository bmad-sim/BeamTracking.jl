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