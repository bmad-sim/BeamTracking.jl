#---------------------------------------------------------------------------------------------------

#@inline 
function RFcavity(tm::SaganCavity, bunch, bmultipoleP, rfP, beamlineP, L)
  species = bunch.species
  p1_over_q_ref = beamlineP.beamline.p_over_q_ref
  rf_omega = rf_omega_calc(rfP, beamlineP.beamline.line[end].s_downstream, species, p1_over_q_ref)
  E1_ref = R_to_E(species, p1_over_q_ref)
  dE_ref = beamlineP.dE_ref
  E0_ref = E1_ref - dE_ref
  p0_over_q_ref = E_to_R(species, E0_ref)
  t_phi0 = rf_phi0_calc(rfP, beamlineP.beamline.species_ref) / rf_omega

  mass = massof(species)
  q = chargeof(species)
  n_cell, L_active = rf_step_calc(tm.n_cell, tm.L_active, rf_omega, L)
  t_ref = 0    # bunch.t_ref ## Relative time tracking assumed for now.
  a = gyromagnetic_anomaly(species)

  if L_active == 0
    if isnothing(bmultipoleP)
      return KernelCall(BeamTracking.sagan_cavity_thin!, (tm.radiation_damping_on, 
                         mass, q, E0_ref, dE_ref,
                         t_ref, emptySA, emptySA, emptySA, a, q*rfP.voltage, rf_omega, t_phi0, L))
    else
      BnL, BsL = get_integrated_strengths(bmultipoleP, L, p0_over_q_ref)
      return KernelCall(BeamTracking.sagan_cavity_thin!, (tm.radiation_damping_on, 
                         mass, q, E0_ref, dE_ref,
                         t_ref, bmultipoleP.order, BnL, BsL, a, q*rfP.voltage, rf_omega, t_phi0, L))
    end

  else
    if isnothing(bmultipoleP)
      return KernelCall(BeamTracking.sagan_cavity_thick!, (tm.radiation_damping_on, rfP.traveling_wave, 
                mass, q, E0_ref, dE_ref, t_ref, n_cell, emptySA, emptySA, emptySA, a,
                q*rfP.voltage, rf_omega, t_phi0, L_active, L))
    else
      Bn, Bs = get_strengths(bmultipoleP, L, p0_over_q_ref)
      return KernelCall(BeamTracking.sagan_cavity_thick!, (tm.radiation_damping_on, rfP.traveling_wave, 
                mass, q, E0_ref, dE_ref, t_ref, n_cell, bmultipoleP.order, Bn, Bs, a,
                q*rfP.voltage, rf_omega, t_phi0, L_active, L))
    end
  end
end
