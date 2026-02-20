#---------------------------------------------------------------------------------------------------

@inline function RFcavity(tm::SaganCavity, bunch, bmultipoleP, rfP, beamlineP, L)
  species = bunch.species
  mass = massof(species)
  q = chargeof(species)

  p1_over_q_ref = beamlineP.beamline.p_over_q_ref
  rf_omega = rf_omega_calc(rfP, beamlineP.beamline.line[end].s_downstream, species, p1_over_q_ref)
  E1_ref = R_to_E(species, p1_over_q_ref)
  dE_ref = beamlineP.dE_ref
  E0_ref = E1_ref - dE_ref
  P0c = sqrt(E0_ref^2 - mass^2)
  p0q = E_to_R(species, E0_ref)
  t_phi0 = rf_phi0_calc(rfP, beamlineP.beamline.species_ref) / rf_omega

  num_cells, L_active = rf_step_calc(tm.num_cells, tm.L_active, rf_omega, L)
  L_active <= L * (1 + eps(L)) || error("Cavity cannot have L_active ($L_active) greater than L ($L)." )
  t_ref = 0    # bunch.t_ref ## Relative time tracking assumed for now.
  a = gyromagnetic_anomaly(species)

  #

  if L == 0
    if isnothing(bmultipoleP)
      return KernelCall(BeamTracking.sagan_cavity_zero_L!,
                   (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), mass, q, P0c, dE_ref,
                   t_ref, Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0))
    else
      m_order = bmultipoleP.order
      KnL, KsL = get_integrated_strengths(bmultipoleP, L, p0q)
      has_mult = (length(m_order) > 0)
      return KernelCall(BeamTracking.sagan_cavity_zero_L!, 
                   (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), m_order, KnL.*p0q, KsL.*p0q, a, q*rfP.voltage, rf_omega, t_phi0))
    end

  elseif L_active == 0
    if isnothing(bmultipoleP)
      return KernelCall(BeamTracking.sagan_cavity_zero_L_active!,
                   (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), mass, q, P0c, dE_ref,
                   t_ref, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0, L))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return KernelCall(BeamTracking.sagan_cavity_zero_L_active!, 
                   (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a, q*rfP.voltage, rf_omega, t_phi0, L))
    end

  else
    if isnothing(bmultipoleP)
      return KernelCall(BeamTracking.sagan_cavity_thick!, 
                (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, num_cells, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], 
                a, q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return KernelCall(BeamTracking.sagan_cavity_thick!, 
                (Val{tm.radiation_damping_on}(), Val{tm.radiation_fluctuations_on}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, num_cells, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a,
                q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L))
    end
  end
end
