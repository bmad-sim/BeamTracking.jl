#

@inline function aperture(tm, bunch, apertureparams, entering::Bool)
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
