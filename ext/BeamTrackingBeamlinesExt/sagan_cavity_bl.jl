#---------------------------------------------------------------------------------------------------

@inline function sagan_cavity(tm::SaganCavity, bunch, bmultipole, rf, E0_ref, dE_ref, rf_omega, t_phi0, L)
  mass = massof(bunch.species)
  q    = chargeof(bunch.species)
  n_cell, L_active = rf_step_calc(tm.n_cell, tm.L_active, rf_omega, L)
  t_ref = 0    # bunch.t_ref ## Relative time tracking assumed for now.
  p_over_q_ref = E_to_R(bunch.species, E0_ref+dE_ref)
  a = gyromagnetic_anomaly(bunch.species)

  if L_active == 0
    if isnothing(bmultipole)
      return KernelCall(BeamTracking.sagan_cavity_thin!, (tm.radiation_damping_on, 
                         mass, q, E0_ref, dE_ref,
                         t_ref, emptySA, emptySA, emptySA, a, q*rf.voltage, rf_omega, t_phi0, L))
    else
      BnL, BsL = get_integrated_strengths(bmultipole, L, p_over_q_ref)
      print(typeof(BnL))
      return KernelCall(BeamTracking.sagan_cavity_thin!, (tm.radiation_damping_on, 
                         mass, q, E0_ref, dE_ref,
                         t_ref, bmultipole.order, BnL, BsL, a, q*rf.voltage, rf_omega, t_phi0, L))
    end

  else
    if isnothing(bmultipole)
      return KernelCall(BeamTracking.sagan_cavity_thick!, (tm.radiation_damping_on, rf.traveling_wave, 
                mass, q, E0_ref, dE_ref, t_ref, n_cell, emptySA, emptySA, emptySA, a,
                q*rf.voltage, rf_omega, t_phi0, L_active, L))
    else
      Bn, Bs = get_strengths(bmultipole, L, p_over_q_ref)
      return KernelCall(BeamTracking.sagan_cavity_thick!, (tm.radiation_damping_on, rf.traveling_wave, 
                mass, q, E0_ref, dE_ref, t_ref, n_cell, bmultipole.order, Bn, Bs, a,
                q*rf.voltage, rf_omega, t_phi0, L_active, L))
    end
  end
end
