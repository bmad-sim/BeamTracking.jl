#---------------------------------------------------------------------------------------------------
# Here both L_active and L are zero

function sagan_cavity_zero_L!(i, coords::Coords, 
                                      val_rad_damping_on, val_rad_fluctuations_on,
                                      mass, q, P0c, dE_ref, t_ref, 
                                      ::Val{has_mult}, m_order, BnL, BsL, 
                                      a, q_voltage, rf_omega, t_phi0) where {has_mult}
@inbounds begin
  # Multipole kick
  if has_mult
    f = q * C_LIGHT / (2 * P0c)    # q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c)
  end

  # Energy kick
  sagan_cavity_kick!(i, coords, a, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)

  # Reference energy shift
  dP0c = dpc_given_dE(P0c, dE_ref, mass)
  reference_momentum_shift!(i, coords, P0c, dP0c, true)
  P0c += dP0c

  # Multipole kick
  if has_mult
    f = q * C_LIGHT / (2 * P0c)    # q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c)
  end
end
end

#---------------------------------------------------------------------------------------------------
# Here L_active is zero and L is nonzero. Note: The case L_active > L is not allowed.

@makekernel fastgtpsa=true function sagan_cavity_zero_L_active!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on,
                        mass, q, P0c, dE_ref, t_ref, 
                        val_has_mult, val_has_sol, m_order, Bn, Bs, a, q_voltage, rf_omega, t_phi0, L)

  # Outside Drift
  f = q * C_LIGHT / P0c    # q_over_p_ref
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
                    val_has_mult, val_has_sol, m_order, Bn .* f, Bs .* f, q, a, mass, P0c, L/2)

  # Energy kick
  sagan_cavity_kick!(i, coords, a, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)

  # Reference energy shift
  dP0c = dpc_given_dE(P0c, dE_ref, mass)
  reference_momentum_shift!(i, coords, P0c, dP0c, true)
  P0c += dP0c

  # Outside Drift
  f = q * C_LIGHT / P0c    # q_over_p_ref
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
                    val_has_mult, val_has_sol, m_order, Bn .* f, Bs .* f, q, a, mass, P0c, L/2)
end

#---------------------------------------------------------------------------------------------------
# Here both L_active and L are non-zero

@makekernel fastgtpsa=true function sagan_cavity_thick!(i, coords::Coords, 
              val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave, 
              mass, q, P0c, dE_ref, t_ref, n_cell, 
              val_has_mult, val_has_sol, m_order, Bn, Bs, 
              a, q_gradient, rf_omega, t_phi0, L_active, L)

  # Outside Drift
  q_over_p_ref = q * C_LIGHT / P0c
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
          val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
          q, a, mass, P0c, (L-L_active)/2)

  # Fringe kick at beginning. 
  sagan_cavity_fringe!(i, coords, a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, +1)

  # Body Loop
  # n_cell == 0 => single kick in center

  if n_cell == 0
    sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
              val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
              q, q_gradient, a, mass, P0c, L_active/2)
    sagan_cavity_kick!(i, coords, a, q_gradient*L_active, rf_omega, t_phi0, t_ref, mass, P0c)
    dP0c = dpc_given_dE(P0c, dE_ref, mass)
    reference_momentum_shift!(i, coords, P0c, dP0c, true)
    P0c += dP0c
    q_over_p_ref = q * C_LIGHT / P0c

    sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
                val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
                q, q_gradient, a, mass, P0c, L_active/2)

  else
    for i_step = 0:n_cell
      i_step == 0 || i_step == n_cell ? kick_factor = 2 : kick_factor = 1

      # Longitudinal kick
      sagan_cavity_kick!(i, coords, a, q_gradient*L_active/(n_cell*kick_factor), rf_omega, t_phi0, t_ref, mass, P0c)

      # Reference energy shift
      dP0c = dpc_given_dE(P0c, dE_ref/(n_cell*kick_factor), mass)
      reference_momentum_shift!(i, coords, P0c, dP0c, true)
      P0c += dP0c
      q_over_p_ref = q * C_LIGHT / P0c

      # Drift
      if i_step == n_cell; break; end
      sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
           val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
           q, q_gradient, a, mass, P0c, L_active/n_cell)
    end
  end

  # Fringe kick at end
  sagan_cavity_fringe!(i, coords, a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, -1)

  # Outside Drift
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
              val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
              q, a, mass, P0c, (L-L_active)/2)
end

#---------------------------------------------------------------------------------------------------
# The "inside" drift is the same as the "outside" drift with the addition of the pondermotive kick.

function sagan_cavity_inside_drift!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on, ::Val{traveling_wave}, 
                        val_has_mult, val_has_sol, m_order, Kn, Ks, 
                        q, q_gradient, a, mass, P0c, L) where {traveling_wave}
  # 1/2 pondermotive kick
  if !traveling_wave
    rf_pondermotive_kick!(i, coords, q_gradient, a, mass, P0c, L/2)
  end

  # Drift
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                    val_has_mult, val_has_sol, m_order, Kn, Ks, q, a, mass, P0c, L)

  # 1/2 pondermotive kick
  # The pondermotive force only occurs if there is a EM wave in the opposite direction from the direction of travel.
  if !traveling_wave
    rf_pondermotive_kick!(i, coords, q_gradient, a, mass, P0c, L/2)
  end
end

#---------------------------------------------------------------------------------------------------
# The "outside" drift does not have the pondermotive kick that the "inside" drift does.

function sagan_cavity_outside_drift!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on, 
                        ::Val{has_mult}, ::Val{has_sol}, m_order, Kn, Ks, 
                        q, a, mass, P0c, L) where {has_mult, has_sol}
  beta0 = P0c / sqrt(P0c^2 + mass^2)
  gamma0 = P0c / (mass * beta0)

  # 1/2 Multipole kick
  if has_mult
    multipole_kick_with_rad!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L/2)
  end

  # Solenoid or Drift
  if has_sol
    exact_solenoid!(i, coords, Kn[1], beta0, gamma0*gamma0, mass/P0c, L)
  else
    exact_drift!(i, coords, beta0, gamma0*gamma0, mass/P0c, L)
  end

  # 1/2 Multipole kick
  if has_mult
    multipole_kick_with_rad!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L/2)
  end
end

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function sagan_cavity_kick!(i, coords::Coords, 
                            a, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  pz = v[i,PZI]
  Pc = (1 + pz) * P0c
  m_over_pc = mass / ((1 + v[i,PZI]) * P0c)
  beta = 1 / sqrt(1 + m_over_pc * m_over_pc)
  t = t_phi0 + t_ref - v[i,ZI] / (beta * C_LIGHT)
  dE = q_voltage * cos(rf_omega * t)

  # 

  if !isnothing(coords.q)
    e_field = (0, 0, dE * C_LIGHT / P0c)
    b_field = (0, 0, 0)
    a_potential = 0
    g_bend = 0
    rotate_spin_field!(i, coords, a, g_bend, mass/P0c, a_potential, a_potential, e_field, b_field, 1)
  end

  # Energy kick

  to_energy_coords!(i, coords, beta)
  rad = dE*dE + 2*sqrt(Pc*Pc + mass^2) * dE + Pc*Pc
  coords.state[i] = vifelse(rad < 0, STATE_LOST_PZ, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  sqrt_rad = vifelse(alive, sqrt(abs(rad)), one(rad))
  pz = vifelse(alive, pz + (sqrt_rad - Pc)/P0c, pz)  

  v[i,PZI] = vifelse(alive, v[i,PZI] + dE/P0c, v[i,PZI])
  to_momentum_coords!(i, coords, mass, P0c, pz)

  #

  if !isnothing(coords.q)
    rotate_spin_field!(i, coords, a, g_bend, mass/P0c, a_potential, a_potential, e_field, b_field, 1)
  end
end

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function rf_pondermotive_kick!(i, coords, q_gradient, a, mass, P0c, L)
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
    sagan_cavity_fringe!(i, coords::Coords, a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, edge)

Fringe kick due to forward traveling wave.

edge = +1 => entering, edge = -1 => exiting
"""
@makekernel fastgtpsa=true function sagan_cavity_fringe!(i, coords,
                                      a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, edge)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  pz = v[i,PZI]
  Pc = (1 + pz) * P0c

  m_over_pc = mass / ((1 + v[i,PZI]) * P0c)
  beta = 1 / sqrt(1 + m_over_pc * m_over_pc)
  t = t_phi0 + t_ref - v[i,ZI] / (beta * C_LIGHT)
  phase = rf_omega * t
  sin_phase, cos_phase = sincos(phase)
  ez_field = q_gradient * cos_phase
  dez_dz_field = q_gradient * sin_phase * rf_omega / C_LIGHT
  dE = -edge * dez_dz_field * (v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]) / 4
  rad = dE*dE + 2*sqrt(Pc*Pc + mass^2) * dE + Pc*Pc
  coords.state[i] = vifelse(rad < 0, STATE_LOST_PZ, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  sqrt_rad = vifelse(alive, sqrt(abs(rad)), one(rad))
  pz = vifelse(alive, pz + (sqrt_rad - Pc)/P0c, pz)  
  # Spin

  if !isnothing(coords.q)
    f = -edge * ez_field * C_LIGHT / (4 * P0c)
    e_field = (f*v[i,XI], f*v[i,YI], 0)
    b_field = (0, 0, 0)
    a_potential = 0
    g_bend = 0
    rotate_spin_field!(i, coords, a, g_bend, mass/P0c, a_potential, a_potential, e_field, b_field, 1)
  end

  # Fringe

  to_energy_coords!(i, coords, beta)
  f = edge * ez_field / (2 * P0c)
  v[i,PXI] = vifelse(alive, v[i,PXI] - f * v[i,XI], v[i,PXI])
  v[i,PYI] = vifelse(alive, v[i,PYI] - f * v[i,YI], v[i,PYI])
  ## Note: v[i,PZI] will be set in to_momentum_coords!
  to_momentum_coords!(i, coords, mass, P0c, pz)

  # Spin
  if !isnothing(coords.q)
    rotate_spin_field!(i, coords, a, g_bend, mass/P0c, a_potential, a_potential, e_field, b_field, 1)
  end
end

#---------------------------------------------------------------------------------------------------


function multipole_kick_with_rad!(i, coords::Coords, 
                  ::Val{rad_damping_on}, ::Val{rad_fluctuations_on}, 
                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L) where {rad_damping_on, rad_fluctuations_on}
  L2 = L / 2
  if rad_damping_on
    deterministic_radiation!(i, coords, q, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
  end

  if rad_fluctuations_on
    E0 = sqrt(P0c^2 + mass^2)
    stochastic_radiation!(i, coords, q, mass, E0, 0, 0, m_order, Kn, Ks, L2)
  end

  multipole_and_spin_kick!(i, coords, m_order, Kn,  Ks, a, mass/P0c, L)

  if rad_fluctuations_on
    stochastic_radiation!(i, coords, q, mass, E0, 0, 0, m_order, Kn, Ks, L2)
  end

  if rad_damping_on
    deterministic_radiation!(i, coords, q, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
  end
end

#---------------------------------------------------------------------------------------------------

"""
    to_energy_coords!(i, coords::Coords, beta)

Convert from `(z, pz)` to `(c(t_ref - t), E/P0c)` coords.
"""
@makekernel fastgtpsa=true function to_energy_coords!(i, coords::Coords, beta)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
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

function dpc_given_dE(old_pc, dE, mass)
  rad = dE*dE + 2*sqrt(old_pc*old_pc + mass^2) * dE
  return rad/(sqrt(rad + old_pc*old_pc) + old_pc)
end
