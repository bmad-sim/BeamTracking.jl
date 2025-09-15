using BeamTracking
using Beamlines

@eles begin
  d1 = Drift(L=2, x_offset = 0.1, tracking_method = Exact())
  d2 = Drift(L=2, tracking_method = Exact())
end

bline_d1 = Beamline([d1], R_ref = 1.0, species_ref=Species("electron"))
bline_d2 = Beamline([d2], R_ref = 1.0, species_ref=Species("electron"))

;
