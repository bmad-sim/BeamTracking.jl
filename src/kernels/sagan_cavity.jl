#---------------------------------------------------------------------------------------------------
# Here both L_active and L are zero

function sagan_cavity_zero_L!(i, coords::Coords, 
                                      val_rad_damping_on, val_rad_fluctuations_on,
                                      mass, q, P0c, dE_ref, t_ref, 
                                      ::Val{has_mult}, m_order, BnL, BsL, 
                                      a, q_voltage, rf_omega, t_phi0, design_params=nothing) where {has_mult}
@inbounds begin
  # Multipole kick
  if has_mult
    f = q * c_light(eltype(coords.v)) / (2 * P0c)    # q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c)
  end

  # Energy kick
  sagan_cavity_kick!(i, coords, a, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)

  # Reference energy shift
  dP0c = dpc_given_dE(P0c, dE_ref, mass)
  reference_momentum_shift!(i, coords, P0c, dP0c, Val{true}())
  P0c += dP0c

  # Multipole kick
  if has_mult
    f = q * c_light(eltype(coords.v)) / (2 * P0c)    # q_over_p_ref / 2
    multipole_and_spin_kick!(i, coords, m_order, BnL .* f, BsL .* f, a, mass/P0c)
  end
end
end

#---------------------------------------------------------------------------------------------------
# A slice of a SaganCavity (n_slices > 1) tracks the part of the one cavity between two cuts: the 
# kicks, fringes and drifts of the whole cavity stay where they are, and a drift that crosses a cut 
# is split. `window` (see sagan_window in the Beamlines extension) says which parts are tracked:
#
#   ka:kb                       the kicks in the window (0-based; empty if kb < ka). The reference 
#                               energy steps of the kicks before it are applied to P0c (and P0c_d) only.
#   outer_in, outer_out         whether the drifts outside the active region are (partly) in the window;
#   pre, post                   whether the inside drift before kick ka and after kick kb are;
#   fringe_in, fringe_out       whether the fringe kicks are;
#   *_full                      whether that drift is entirely in the window: its nominal length is then 
#                               computed here (bitwise as without slices), otherwise *_len is its length;
#   *_kick_in, *_kick_out       whether the half kicks at the start and end of that drift are in the window 
#                               (they keep the strength of the whole drift, see sagan_cavity_outside_drift!).
#
# For n_slices = 1 the window is the whole cavity. For a Bunch species other than the reference species 
# (design_params = (mass, P0c, dE_ref) of the reference species, else nothing), the RF phase is relative 
# to the design particle: t_ref includes the lag of the Bunch reference particle behind it since the 
# start of the window (the time shift at the end of each slice re-references z to the design particle).

@inline sagan_lag(t_ref, L, P0c, mass, P0c_d, mass_d, c) = t_ref + L*(sqrt(P0c^2 + mass^2)/P0c - sqrt(P0c_d^2 + mass_d^2)/P0c_d)/c

#---------------------------------------------------------------------------------------------------
# Here L_active is zero and L is nonzero: one kick at L/2. Note: The case L_active > L is not allowed.

@makekernel fastgtpsa=true function sagan_cavity_zero_L_active!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on,
                        mass, q, P0c, dE_ref, t_ref, 
                        val_has_mult, val_has_sol, m_order, Bn, Bs, a, q_voltage, rf_omega, t_phi0, L, window, design_params=nothing)
  c = c_light(eltype(coords.v))
  if !isnothing(design_params)
    mass_d, P0c_d, dE_d = design_params
  end

  # The kick is before the window
  if window.ka > 0
    P0c += dpc_given_dE(P0c, dE_ref, mass)
  end

  # Outside Drift
  if window.outer_in
    d = window.outer_in_full ? L/2 : window.outer_in_len
    f = q * c / P0c    # q_over_p_ref
    sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
                      val_has_mult, val_has_sol, m_order, Bn .* f, Bs .* f, q, a, mass, P0c, d, 
                      window.outer_in_kick_in, window.outer_in_kick_out, L/2)
    if !isnothing(design_params)
      t_ref = sagan_lag(t_ref, d, P0c, mass, P0c_d, mass_d, c)
    end
  end

  if window.kb >= window.ka
    # Energy kick
    sagan_cavity_kick!(i, coords, a, q_voltage, rf_omega, t_phi0, t_ref, mass, P0c)

    # Reference energy shift
    dP0c = dpc_given_dE(P0c, dE_ref, mass)
    reference_momentum_shift!(i, coords, P0c, dP0c, Val{true}())
    P0c += dP0c
  end

  # Outside Drift
  if window.outer_out
    d = window.outer_out_full ? L/2 : window.outer_out_len
    f = q * c / P0c    # q_over_p_ref
    sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
                      val_has_mult, val_has_sol, m_order, Bn .* f, Bs .* f, q, a, mass, P0c, d, 
                      window.outer_out_kick_in, window.outer_out_kick_out, L/2)
  end
end

#---------------------------------------------------------------------------------------------------
# Here both L_active and L are non-zero. n_cells == 0 => single kick in the center of the active region, 
# else n_cells + 1 kicks (half kicks at the ends of the active region).

@makekernel fastgtpsa=true function sagan_cavity_thick!(i, coords::Coords, 
              val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave, 
              mass, q, P0c, dE_ref, t_ref, n_cells, 
              val_has_mult, val_has_sol, m_order, Bn, Bs, 
              a, q_gradient, rf_omega, t_phi0, L_active, L, window, design_params=nothing)
  c = c_light(eltype(coords.v))
  if !isnothing(design_params)
    mass_d, P0c_d, dE_d = design_params
  end
  # Length of the inside drifts (between kicks; on either side of the kick for n_cells == 0)
  L_cell = n_cells == 0 ? L_active/2 : L_active/n_cells

  # The reference energy steps of the kicks before the window
  for i_step = 0:window.ka-1
    kick_div = n_cells == 0 ? 1 : (i_step == 0 || i_step == n_cells ? 2*n_cells : n_cells)
    P0c += dpc_given_dE(P0c, dE_ref/kick_div, mass)
    if !isnothing(design_params)
      P0c_d += dpc_given_dE(P0c_d, dE_d/kick_div, mass_d)
    end
  end
  q_over_p_ref = q * c / P0c

  # Outside Drift
  if window.outer_in
    d = window.outer_in_full ? (L-L_active)/2 : window.outer_in_len
    sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
            val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
            q, a, mass, P0c, d, window.outer_in_kick_in, window.outer_in_kick_out, (L-L_active)/2)
    if !isnothing(design_params)
      t_ref = sagan_lag(t_ref, d, P0c, mass, P0c_d, mass_d, c)
    end
  end

  # Fringe kick at beginning. 
  if window.fringe_in
    sagan_cavity_fringe!(i, coords, a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, +1)
  end

  # Body: the inside drift before the first kick of the window, the kicks with the drifts between them, 
  # and the inside drift after the last kick
  if window.pre
    d = window.pre_full ? L_cell : window.pre_len
    sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
         val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
         q, q_gradient, a, mass, P0c, d, window.pre_kick_in, window.pre_kick_out, L_cell)
    if !isnothing(design_params)
      t_ref = sagan_lag(t_ref, d, P0c, mass, P0c_d, mass_d, c)
    end
  end

  for i_step = window.ka:window.kb
    kick_div = n_cells == 0 ? 1 : (i_step == 0 || i_step == n_cells ? 2*n_cells : n_cells)

    # Longitudinal kick
    sagan_cavity_kick!(i, coords, a, q_gradient*L_active/kick_div, rf_omega, t_phi0, t_ref, mass, P0c)

    # Reference energy shift
    dP0c = dpc_given_dE(P0c, dE_ref/kick_div, mass)
    reference_momentum_shift!(i, coords, P0c, dP0c, Val{true}())
    P0c += dP0c
    q_over_p_ref = q * c / P0c
    if !isnothing(design_params)
      P0c_d += dpc_given_dE(P0c_d, dE_d/kick_div, mass_d)
    end

    # Drift
    if i_step == window.kb; break; end
    sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
         val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
         q, q_gradient, a, mass, P0c, L_cell)
    if !isnothing(design_params)
      t_ref = sagan_lag(t_ref, L_cell, P0c, mass, P0c_d, mass_d, c)
    end
  end

  if window.post
    d = window.post_full ? L_cell : window.post_len
    sagan_cavity_inside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, val_traveling_wave,
         val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
         q, q_gradient, a, mass, P0c, d, window.post_kick_in, window.post_kick_out, L_cell)
    if !isnothing(design_params)
      t_ref = sagan_lag(t_ref, d, P0c, mass, P0c_d, mass_d, c)
    end
  end

  # Fringe kick at end
  if window.fringe_out
    sagan_cavity_fringe!(i, coords, a, q_gradient, rf_omega, t_phi0, t_ref, mass, P0c, -1)
  end

  # Outside Drift
  if window.outer_out
    d = window.outer_out_full ? (L-L_active)/2 : window.outer_out_len
    sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on,
                val_has_mult, val_has_sol, m_order, Bn .* q_over_p_ref, Bs .* q_over_p_ref, 
                q, a, mass, P0c, d, window.outer_out_kick_in, window.outer_out_kick_out, (L-L_active)/2)
  end
end

#---------------------------------------------------------------------------------------------------
# The "inside" drift is the same as the "outside" drift with the addition of the pondermotive kick.

function sagan_cavity_inside_drift!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on, ::Val{traveling_wave}, 
                        val_has_mult, val_has_sol, m_order, Kn, Ks, 
                        q, q_gradient, a, mass, P0c, L, kick_in=true, kick_out=true, L_kick=L) where {traveling_wave}
  # 1/2 pondermotive kick
  if !traveling_wave && kick_in
    rf_pondermotive_kick!(i, coords, q_gradient, a, mass, P0c, L_kick/2)
  end

  # Drift
  sagan_cavity_outside_drift!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                    val_has_mult, val_has_sol, m_order, Kn, Ks, q, a, mass, P0c, L, kick_in, kick_out, L_kick)

  # 1/2 pondermotive kick
  # The pondermotive force only occurs if there is a EM wave in the opposite direction from the direction of travel.
  if !traveling_wave && kick_out
    rf_pondermotive_kick!(i, coords, q_gradient, a, mass, P0c, L_kick/2)
  end
end

#---------------------------------------------------------------------------------------------------
# The "outside" drift does not have the pondermotive kick that the "inside" drift does.
# The drift of length L is a part of a drift of length L_kick (split by a slice cut): the half kicks 
# at its ends keep their strength and are applied by the slice that has that end (kick_in, kick_out).

function sagan_cavity_outside_drift!(i, coords::Coords, 
                        val_rad_damping_on, val_rad_fluctuations_on, 
                        ::Val{has_mult}, ::Val{has_sol}, m_order, Kn, Ks, 
                        q, a, mass, P0c, L, kick_in=true, kick_out=true, L_kick=L) where {has_mult, has_sol}
  beta0 = P0c / sqrt(P0c^2 + mass^2)
  gamma0 = P0c / (mass * beta0)

  # 1/2 Multipole kick
  if has_mult && kick_in
    multipole_kick_with_rad!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L_kick/2)
  end

  # Solenoid or Drift
  if has_sol
    exact_solenoid!(i, coords, Kn[1], beta0, gamma0*gamma0, mass/P0c, nothing, L)
  else
    exact_drift!(i, coords, 0, beta0, gamma0*gamma0, mass/P0c, L)
  end

  # 1/2 Multipole kick
  if has_mult && kick_out
    multipole_kick_with_rad!(i, coords, val_rad_damping_on, val_rad_fluctuations_on, 
                                  q, m_order, Kn, Ks, a, mass, P0c, beta0, L_kick/2)
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
  t = t_phi0 + t_ref - v[i,ZI] / (beta * c_light(eltype(coords.v)))
  dE = q_voltage * cos(rf_omega * t)

  # 

  if !isnothing(coords.q)
    e_field = (0, 0, dE * c_light(eltype(coords.v)) / P0c)
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
  t = t_phi0 + t_ref - v[i,ZI] / (beta * c_light(eltype(coords.v)))
  phase = rf_omega * t
  sin_phase, cos_phase = sincos(phase)
  ez_field = q_gradient * cos_phase
  dez_dz_field = q_gradient * sin_phase * rf_omega / c_light(eltype(coords.v))
  dE = -edge * dez_dz_field * (v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]) / 4
  rad = dE*dE + 2*sqrt(Pc*Pc + mass^2) * dE + Pc*Pc
  coords.state[i] = vifelse(rad < 0, STATE_LOST_PZ, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  sqrt_rad = vifelse(alive, sqrt(abs(rad)), one(rad))
  pz = vifelse(alive, pz + (sqrt_rad - Pc)/P0c, pz)  
  # Spin

  if !isnothing(coords.q)
    f = -edge * ez_field * c_light(eltype(coords.v)) / (4 * P0c)
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
    deterministic_radiation_multipole!(i, coords, q, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
  end

  if rad_fluctuations_on
    E0 = sqrt(P0c^2 + mass^2)
    stochastic_radiation!(i, coords, q, mass, E0, 0, nothing, nothing, m_order, Kn, Ks, L2)
  end

  multipole_and_spin_kick!(i, coords, m_order, Kn,  Ks, a, mass/P0c, L)

  if rad_fluctuations_on
    stochastic_radiation!(i, coords, q, mass, E0, 0, nothing, nothing, m_order, Kn, Ks, L2)
  end

  if rad_damping_on
    deterministic_radiation_multipole!(i, coords, q, mass, P0c/beta0, 0, m_order, Kn, Ks, L2)
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
