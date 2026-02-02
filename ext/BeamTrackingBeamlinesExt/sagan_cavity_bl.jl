#---------------------------------------------------------------------------------------------------

@inline function sagan_cavity(tm::SaganCavity, bunch, bmultipole, rf, E0_ref, dE_ref, rf_omega, t_phi0, L)

  mass = massof(bunch.species)
  q    = chargeof(bunch.species)
  n_cell, L_active = rf_step_calc(tm.n_cell, tm.L_active, rf_omega, L)
  t_ref = 0    # bunch.t_ref ## Relative time tracking assumed for now.
  p_over_q_ref = E_to_R(bunch.species, E0_ref+dE_ref)
  a = gyromagnetic_anomaly(bunch.species)

  if L_active == 0
    m_order, Bnl, Bsl, Bsol = get_multipole_fields(bmultipole, L, p_over_q_ref)
    return KernelCall(BeamTracking.sagan_cavity_thin!, (tm.radiation_damping_on, 
                         mass, q, E0_ref, dE_ref,
                         t_ref, m_order, Bnl, Bsl, Bsol, a, q*rf.voltage, rf_omega, t_phi0, L))

  else
    m_order, Bn, Bs, Bsol = get_multipole_fields(bmultipole, L, p_over_q_ref)
    return KernelCall(BeamTracking.sagan_cavity_thick!, (tm.radiation_damping_on, rf.traveling_wave, 
                mass, q, E0_ref, dE_ref, t_ref, n_cell, m_order, Bn, Bs, Bsol, a,
                q*rf.voltage, rf_omega, t_phi0, L_active, L))
  end
end
