using BeamTracking
using Beamlines

@elements begin
  d_error = LineElement(L=1, x1_limit=1, dx=5)

  d1_rect = Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_active = false)
  d2_rect = Drift(L=1,               x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_at = ApertureAt.Exit)
  d3_rect = Drift(L=1,               x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_shape = ApertureShape.Rectangular, aperture_at = ApertureAt.BothEnds)
  d4_rect = Drift(L=1, x1_limit = 1, x2_limit = 2,                             aperture_shape = ApertureShape.Rectangular)

  d1_ellip =  Drift(L=1, x1_limit = -Inf, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_active = false)
  d2_ellip =  Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5, aperture_at = ApertureAt.Exit)
  d3_ellip =  Drift(L=1, x1_limit = 1, x2_limit = 2, y1_limit = 3, y2_limit = 5)
  d4_ellip =  Drift(L=1, x1_limit = 1, x2_limit = Inf, y1_limit = 3, y2_limit = 5)
end
b_error = Beamline([d_error], pc_ref = 1e6, species_ref=Species("electron"))

b1r = Beamline([d1_rect], pc_ref = 1e6, species_ref=Species("electron"))
b2r = Beamline([d2_rect], pc_ref = 1e6, species_ref=Species("electron"))
b3r = Beamline([d3_rect], pc_ref = 1e6, species_ref=Species("electron"))
b4r = Beamline([d4_rect], pc_ref = 1e6, species_ref=Species("electron"))

b1e = Beamline([d1_ellip], pc_ref = 1e6, species_ref=Species("electron"))
b2e = Beamline([d2_ellip], pc_ref = 1e6, species_ref=Species("electron"))
b3e = Beamline([d3_ellip], pc_ref = 1e6, species_ref=Species("electron"))
b4e = Beamline([d4_ellip], pc_ref = 1e6, species_ref=Species("electron"))

;