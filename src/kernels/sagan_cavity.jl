#---------------------------------------------------------------------------------------------------
# "thin" means L_active is zero

@makekernel fastgtpsa=true function sagan_cavity_thin!(i, coords::Coords, radiation_damping_on, 
                                      mass, q, E0_ref, dE_ref, t_ref, 
                                      m_order, BnL, BsL, a, q_voltage, rf_omega, t_phi0, L)
  L_out = L / 2           # Length outside of active region
  P0c = sqrt(E0_ref^2 - mass^2)
  q_over_p_ref = q * C_LIGHT / P0c

  # Outside Drift
  if L == 0
    f = q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c, L)
  else
    f = q_over_p_ref / L_out
    sagan_cavity_outside_drift!(i, coords, radiation_damping_on, q,
                                            m_order, BnL .* f, BsL .* f, a, mass, P0c, L_out)
  end

  # Energy kick
  sagan_cavity_kick!(i, coords, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)

  # Reference energy shift
  dP0c = dpc_given_dE(P0c, dE_ref, mass)
  reference_energy_shift!(i, coords, P0c, dP0c)
  P0c += dP0c
  q_over_p_ref = q * C_LIGHT / P0c

  # Outside Drift
  if L == 0
    f = q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c, L)
  else
    f = q_over_p_ref / L_out
    sagan_cavity_outside_drift!(i, coords, radiation_damping_on,  q,
                                            m_order, BnL .* f ,  BsL .* f, a, mass, P0c, L_out)
  end
end

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function sagan_cavity_thick!(i, coords::Coords, 
              radiation_damping_on, traveling_wave, 
              mass, q, E0_ref, dE_ref, t_ref, n_cell, m_order, Bn, Bs, 
              a, q_voltage, rf_omega, t_phi0, L_active, L)

  L_out = (L - L_active) / 2           # Length outside of active region
  q_gradient = q_voltage / L_active    # Effective gradient
  P0c = sqrt(E0_ref^2 - mass^2)
  q_over_p_ref = q * C_LIGHT / P0c

  # Outside Drift
  sagan_cavity_outside_drift!(i, coords, radiation_damping_on, q,
                         m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, a, mass, P0c, L_out)

  # Fringe kick at beginning. 
  sagan_cavity_fringe!(i, coords, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, +1)

  # Body Loop
  # n_cell == 0 => single kick in center

  if n_cell == 0
    sagan_cavity_inside_drift!(i, coords, radiation_damping_on, traveling_wave, q, q_gradient, 
              m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, a, mass, P0c, L_active/2)
    sagan_cavity_kick!(i, coords, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)
    dP0c = dpc_given_dE(P0c, dE_ref, mass)
    reference_energy_shift!(i, coords, P0c, dP0c)
    P0c += dP0c
    q_over_p_ref = q * C_LIGHT / P0c

    sagan_cavity_inside_drift!(i, coords, radiation_damping_on, traveling_wave, q, q_gradient, 
                             m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, a, mass, P0c, L_active/2)

  else
    for i_step = 0:n_cell
      i_step == 0 || i_step == n_cell ? kick_factor = 2 : kick_factor = 1

      # Longitudinal kick
      sagan_cavity_kick!(i, coords, q_voltage/(n_cell*kick_factor), rf_omega, t_phi0, t_ref, mass, P0c)

      # Reference energy shift
      dP0c = dpc_given_dE(P0c, dE_ref/(n_cell*kick_factor), mass)
      reference_energy_shift!(i, coords, P0c, dP0c)
      P0c += dP0c

      # Drift
      if i_step == n_cell; break; end
      sagan_cavity_inside_drift!(i, coords, radiation_damping_on, traveling_wave, q, q_gradient, 
          m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, a, mass, P0c, L_active/n_cell)
    end
  end

  # Fringe kick at end
  sagan_cavity_fringe!(i, coords, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, -1)

  # Outside Drift
  sagan_cavity_outside_drift!(i, coords, radiation_damping_on, q,
                          m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, a, mass, P0c, L_out)
end

#---------------------------------------------------------------------------------------------------
# The "inside" drift is the same as the "outside" drift with the addition of the pondermotive kick.

@makekernel fastgtpsa=true function sagan_cavity_inside_drift!(i, coords::Coords, 
        radiation_damping_on, traveling_wave, q, q_gradient, m_order, Kn, Ks, a, mass, P0c, L)

  beta0 = P0c / sqrt(P0c^2 + mass^2)
  gamma0 = P0c / (mass * beta0)

  # 1/2 pondermotive kick
  if !traveling_wave
    rf_pondermotive_kick!(i, coords, q_gradient, P0c, L/2)
  end

  # Drift
  sagan_cavity_outside_drift!(i, coords, radiation_damping_on, q, 
                                    m_order, Kn, Ks, a, mass, P0c, L)

  # 1/2 pondermotive kick
  # The pondermotive force only occurs if there is a EM wave in the opposite direction from the direction of travel.
  if !traveling_wave
    rf_pondermotive_kick!(i, coords, q_gradient, P0c, L/2)
  end
end

#---------------------------------------------------------------------------------------------------
# The "outside" drift does not have the pondermotive kick that the "inside" drift does.

@makekernel fastgtpsa=true function sagan_cavity_outside_drift!(i, coords::Coords, 
                                     radiation_damping_on, q, m_order, Kn, Ks, a, mass, P0c, L)
  beta0 = P0c / sqrt(P0c^2 + mass^2)
  gamma0 = P0c / (mass * beta0)
  length(m_order) > 0 ? has_multipoles = true : has_multipoles = false
  has_sol = (has_multipoles && m_order[1] == 0)

  # 1/2 Multipole kick
  if has_multipoles
    multipole_kick_with_rad!(i, coords, radiation_damping_on,
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L/2)
  end

  # Drift
  if has_sol
    exact_solenoid!(i, coords, Kn[1], beta0, gamma0*gamma0, mass/P0c, L)
  else
    exact_drift!(i, coords, beta0, gamma0*gamma0, mass/P0c, L)
  end

  # 1/2 Multipole kick
  if has_multipoles
    multipole_kick_with_rad!(i, coords, radiation_damping_on,
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L/2)
  end
end

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function sagan_cavity_kick!(i, coords::Coords, 
                            q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  # Energy kick
  pz = v[i,PZI]
  Pc = (1 + pz) * P0c

  to_energy_coords!(i, coords, mass, P0c)

  t = t_phi0 + t_ref - v[i,ZI] / C_LIGHT
  dE = q_voltage * cos(rf_omega * t)

  rad = dE*dE + 2*sqrt(Pc*Pc + mass^2) * dE + Pc*Pc
  coords.state[i] = vifelse(rad < 0, STATE_LOST_PZ, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  sqrt_rad = vifelse(alive, sqrt(abs(rad)), 1.0)
  pz = vifelse(alive, pz + (sqrt_rad - Pc)/P0c, pz)  

  v[i,PZI] = vifelse(alive, v[i,PZI] + dE/P0c, v[i,PZI])
  to_momentum_coords!(i, coords, mass, P0c, pz)
end

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function rf_pondermotive_kick!(i, coords, q_gradient, P0c, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  inv_rel_p = 1 / (1 + v[i,PZI])
  coef = q_gradient^2 * L / (8 * P0c^2) * inv_rel_p

  v[i,PXI] = vifelse(alive, v[i,PXI] - coef * v[i,XI], v[i,PXI])
  v[i,PYI] = vifelse(alive, v[i,PYI] - coef * v[i,YI], v[i,PYI])
  v[i,ZI]  = vifelse(alive, v[i,ZI]  - coef * inv_rel_p * (v[i,XI]*v[i,XI]+v[i,YI]*v[i,YI])/2, v[i,ZI])
end

#---------------------------------------------------------------------------------------------------

"""
    sagan_cavity_fringe!(i, coords::Coords, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, edge)

Fringe kick due to forward traveling wave.

edge = +1 => entering, edge = -1 => exiting
"""
@makekernel fastgtpsa=true function sagan_cavity_fringe!(i, coords,
                                      q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, edge)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  pz = v[i,PZI]
  Pc = (1 + pz) * P0c
  
  to_energy_coords!(i, coords, mass, P0c)
  t = t_phi0 + t_ref - v[i,ZI] / C_LIGHT
  phase = rf_omega * t
  ez_field = q_gradient * cos(phase) 
  dez_dz_field = q_gradient * sin(phase) * rf_omega / C_LIGHT
  dE = -edge * dez_dz_field * (v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]) / 4

  rad = dE*dE + 2*sqrt(Pc*Pc + mass^2) * dE + Pc*Pc
  coords.state[i] = vifelse(rad < 0, STATE_LOST_PZ, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  sqrt_rad = vifelse(alive, sqrt(abs(rad)), 1.0)
  pz = vifelse(alive, pz + (sqrt_rad - Pc)/P0c, pz)  

  f = edge * ez_field / (2 * P0c)
  v[i,PXI] = vifelse(alive, v[i,PXI] - f * v[i,XI], v[i,PXI])
  v[i,PYI] = vifelse(alive, v[i,PYI] - f * v[i,YI], v[i,PYI])
  ## Note: v[i,PZI] will be set in to_momentum_coords!

  to_momentum_coords!(i, coords, mass, P0c, pz)
end

#---------------------------------------------------------------------------------------------------


@makekernel fastgtpsa=true function multipole_kick_with_rad!(i, coords::Coords, 
                  radiation_damping_on, q_charge, m_order, Kn, Ks, a, mass, P0c, beta0, L)
  L2 = L / 2
  if radiation_damping_on
    deterministic_radiation!(i, coords, q_charge, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
  end

  multipole_and_spin_kick!(i, coords, m_order, Kn.*L,  Ks.*L, a, mass/P0c, L)

  if radiation_damping_on
    deterministic_radiation!(i, coords, q_charge, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
  end
end

#---------------------------------------------------------------------------------------------------

"""
    to_energy_coords!(i, coords::Coords, mass, P0c)

Convert from `(z, pz)` to `(c(t_ref - t), E/P0c)` coords.
"""
@makekernel fastgtpsa=true function to_energy_coords!(i, coords::Coords, mass, P0c)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  m_over_pc = mass / ((1 + v[i,PZI]) * P0c)
  beta = 1 / sqrt(1 + m_over_pc * m_over_pc)
  v[i,ZI]  = vifelse(alive, v[i,ZI]/ beta, v[i,ZI])
  v[i,PZI] = vifelse(alive, (1 + v[i,PZI]) / beta, v[i,PZI])
end

#---------------------------------------------------------------------------------------------------

"""
    to_momentum_coords!(i, coords::Coords, mass, P0c, pz)

Convert from `(c(t_ref - t), E/P0c)` to `(z, pz)` coords. `pz` is an input since it can be computed
more accurately by the calling routine.
"""
@makekernel fastgtpsa=true function to_momentum_coords!(i, coords::Coords, mass, P0c, pz)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  m_over_pc = mass / ((1 + pz) * P0c)
  beta = 1 / sqrt(1 + m_over_pc * m_over_pc)
  v[i,ZI]  = vifelse(alive, beta * v[i,ZI], v[i,ZI])
  v[i,PZI] = vifelse(alive, pz, v[i,PZI])
end

#---------------------------------------------------------------------------------------------------
"""
    reference_energy_shift!(i, coords::Coords, P0c_old, dP0c)

Shift coordinates due to a change in reference energy `dE`.
"""
@makekernel fastgtpsa=true function reference_energy_shift!(i, coords::Coords, P0c_old, dP0c)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  P0c_new = P0c_old + dP0c
  P0c_ratio = P0c_old / P0c_new

  v[i,PXI] = vifelse(alive, P0c_ratio * v[i,PXI], v[i,PXI])
  v[i,PYI] = vifelse(alive, P0c_ratio * v[i,PYI], v[i,PYI])
  v[i,PZI] = vifelse(alive, P0c_ratio * v[i,PZI] - dP0c / P0c_new, v[i,PZI])
end

#---------------------------------------------------------------------------------------------------

function dpc_given_dE(old_pc, dE, mass)
  rad = dE*dE + 2*sqrt(old_pc*old_pc + mass^2) * dE + old_pc*old_pc
  if rad < 0; error("Change in reference energy too negative."); end
  return sqrt(rad) - old_pc
end
