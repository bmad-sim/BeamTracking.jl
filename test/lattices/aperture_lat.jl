using BeamTracking
using Beamlines

@eles begin
  d1_rect = Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_active = false, tracking_method = Linear())
  d2_rect = Drift(L=1,               x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_at = ApertureAt.Exit, tracking_method = Linear())
  d3_rect = Drift(L=1,               x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_at = ApertureAt.BothEnds, tracking_method = Linear())
  d4_rect = Drift(L=1, x1_limit = 1, x2_limit = 2,                             aperture_shape = ApertureShape.Rectangular, tracking_method = Linear())

  d1_ellip =  Drift(L=1, x1_limit = -Inf, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_active = false, tracking_method = Linear())
  d2_ellip =  Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_at = ApertureAt.Exit, tracking_method = Linear())
  d3_ellip =  Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5, tracking_method = Linear())
  d4_ellip =  Drift(L=1, x1_limit = 1, x2_limit = Inf, y1_limit = 3, y2_limit = 5, tracking_method = Linear())
end

b1r = Beamline([d1_rect], R_ref = 1.0, species_ref=Species("electron"))
b2r = Beamline([d2_rect], R_ref = 1.0, species_ref=Species("electron"))
b3r = Beamline([d3_rect], R_ref = 1.0, species_ref=Species("electron"))
b4r = Beamline([d4_rect], R_ref = 1.0, species_ref=Species("electron"))

b1e = Beamline([d1_ellip], R_ref = 1.0, species_ref=Species("electron"))
b2e = Beamline([d2_ellip], R_ref = 1.0, species_ref=Species("electron"))
b3e = Beamline([d3_ellip], R_ref = 1.0, species_ref=Species("electron"))
b4e = Beamline([d4_ellip], R_ref = 1.0, species_ref=Species("electron"))

;