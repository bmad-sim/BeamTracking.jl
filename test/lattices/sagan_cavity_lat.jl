using BeamTracking
using Beamlines

@eles begin
  sc1 = RFCavity(L = 2.0, voltage = 0e6, rf_frequency =1e9, dE_ref = 0e6,
                  tracking_method = SaganCavity(n_cell = 2))
  sc2 = RFCavity(L = 2.0, voltage = 0e6, rf_frequency =1e9, dE_ref = 1e6,
                  tracking_method = SaganCavity(n_cell = 2))
  sc3 = RFCavity(L = 2.0, voltage = 0e6, rf_frequency =1e9, dE_ref = 0e6,
                  tracking_method = SaganCavity(n_cell = 2))
  m = Marker(E_ref = 1e7, species_ref = Species("electron"))
end

lat = Lattice([m, sc1, sc2])

;
