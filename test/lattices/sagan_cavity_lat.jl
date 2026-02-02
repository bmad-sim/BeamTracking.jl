using BeamTracking
using Beamlines

species = Species("#3He--")
E0 = 1e10

@elements begin
  sc1 = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.1*E0,
                  tracking_method = SaganCavity(n_cell = 2), traveling_wave = false, phi0 = 0.1)
  sc2 = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.1*E0,
                  tracking_method = SaganCavity(n_cell = 3, L_active = 0.0), zero_phase = PhaseReference.BelowTransition)
  sc3 = RFCavity(L = 0.0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.1*E0,
                  tracking_method = SaganCavity(n_cell = 4, L_active = 0.0), zero_phase = PhaseReference.AboveTransition)
  m = Marker(E_ref = E0, species_ref = species)
end

@elements begin
  sc4 = RFCavity(L = 0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.1*E0, Ksol = 0.01,
                  tracking_method = SaganCavity(n_cell = 2), traveling_wave = false, phi0 = 0.1)
  sc5 = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.1*E0, Ks0 = 0.001,
                  tracking_method = SaganCavity(n_cell = 3), zero_phase = PhaseReference.BelowTransition)
  sc6 = RFCavity(L = 2.0, voltage = 0.2*E0, rf_frequency =1e9, dE_ref = 0.0*E0, 
                  tracking_method = SaganCavity(n_cell = 4), zero_phase = PhaseReference.AboveTransition)
  m2 = Marker(E_ref = E0, species_ref = species)
end

lat = Lattice([m, sc1, sc2, sc3])
lat2 = Lattice([m2, sc4, sc5, sc6])

;
