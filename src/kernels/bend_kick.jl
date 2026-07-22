"""
bkb_multipole!()

This integrator uses Bend-Kick-Bend to track a beam through
a curved, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the quadrupole component.

Arguments
—————————
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'k0'       -- dipole strength
- 'mm'       -- order of multipoles
- 'kn'       -- normal multipole strengths 
- 'ks'       -- skew multipole strengths 
- 'L'        -- length
"""
@makekernel fastgtpsa=true function bkb_multipole!(i, coords::Coords, s, radiation_params, tilde_m, beta_0, a, g, w, w_inv, k0, mm, kn, ks, L)
  knl = kn .* L ./ 2
  ksl = ks .* L ./ 2

  rel_p = 1 + coords.v[i,PZI]
  px = coords.v[i,PXI]
  py = coords.v[i,PYI]
  P_s2 = rel_p*rel_p - px*px - py*py
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, g, tilde_m, mm, kn, ks, 1, L / 2)
  end

  if !isnothing(radiation_params)
    q, mc2, E_ref = radiation_params
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, g, mm, kn, ks, L / 2)
  end

  multipole_kick!(i, coords, mm, knl, ksl, 1)
  exact_bend!(i, coords, g*L, g, k0, tilde_m, beta_0, a, L)
  multipole_kick!(i, coords, mm, knl, ksl, 1)

  if !isnothing(radiation_params)
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, g, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, g, tilde_m, mm, kn, ks, 1, L / 2)
  end

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end
end 

"""
    exact_bend!(i, coords::Coords, theta, g, Kn0, tilde_m, beta_0, a, L)

Tracks a particle through a sector bend via exact tracking.

#Arguments
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'Kn0'      -- normalized dipole field
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'a'        -- anomalous gyromagnetic ratio
- 'L'        -- length
"""
@makekernel fastgtpsa=true function exact_bend!(i, coords::Coords, theta, g, Kn0, tilde_m, beta_0, a, L)
  v = coords.v
  px_0 = v[i,PXI]
  py_0 = v[i,PYI]
  rel_p = 1 + v[i,PZI]
  rel_p2 = rel_p*rel_p

  pt2 = rel_p2 - v[i,PYI]*v[i,PYI]
  good_momenta = (pt2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pt2_1 = one(pt2)
  pt = sqrt(vifelse(good_momenta, pt2, pt2_1))

  arg = v[i,PXI] / pt
  abs_arg = abs(arg)
  arg_1 = one(arg)
  arg_0 = zero(arg)
  good_arg = (abs_arg < arg_1)
  coords.state[i] = vifelse(!good_arg & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)

  phi1 = theta + asin(vifelse(good_arg, arg, arg_0))
  gp = Kn0 / pt
  h = 1 + g*v[i,XI] 
  splus, cplus = sincos(phi1)
  sinc_theta = sincu(theta)
  cosc_theta = one_cos_norm(theta)
  sgn = sign(L)
  alpha_helper = h*L*sinc_theta
  alpha = alpha_helper*(2*splus - gp*alpha_helper)

  cond = cplus*cplus + gp*alpha
  good_cond = (cond > 0)
  coords.state[i] = vifelse(!good_cond & alive, STATE_LOST, coords.state[i]) # particle does not intersect the exit face
  alive = (coords.state[i] == STATE_ALIVE)
  cond_1 = one(cond)
  nasty_sqrt = sqrt(vifelse(good_cond, cond, cond_1))

  gp_0 = zero(gp)
  abs_gp = abs(gp)
  good_gp = (abs_gp > gp_0)
  gp_1 = one(gp)
  gp_safe = vifelse(good_gp, gp, gp_1)
  pos_cplus = (cplus > 0)
  xi1 = alpha/(nasty_sqrt + cplus)
  xi2 = (nasty_sqrt - cplus)/gp_safe
  xi = vifelse(!good_gp | pos_cplus, xi1, xi2)

  Lcv = -sgn*(L*sinc_theta + v[i,XI]*sin(theta)) 
  negative_Lcv = -Lcv
  thetap = 2*(phi1 - sgn*atan2(xi, negative_Lcv)) 
  Lp = sgn*sqrt(Lcv*Lcv + xi*xi)/sincu(thetap/2) 

  new_x = v[i,XI]*cos(theta) - L*L*g*cosc_theta + xi
  new_px = pt*sin(phi1 - thetap)
  new_y = v[i,YI] + v[i,PYI]*Lp/pt
  new_z = v[i,ZI] - rel_p*Lp/pt + L*rel_p/sqrt(tilde_m*tilde_m + rel_p*rel_p)/beta_0
  v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    ps2 = pt2 - px_0*px_0
    good_momenta = (ps2 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

    beta_gamma = rel_p/tilde_m
    beta_gamma2 = beta_gamma*beta_gamma
    gamma_minus_1 = beta_gamma2/(1 + sqrt(1 + beta_gamma2))
    gamma = gamma_minus_1 + 1
    vt_over_rel_p = Lp/pt
    factor = vt_over_rel_p*Kn0
    coeff = a*factor*gamma_minus_1*py_0/rel_p2

    o1 = coeff*px_0
    o2 = coeff*py_0 - a*gamma*factor
    o3 = coeff*ps
    q1 = expq((o1, o2, o3), alive)

    theta_0 = zero(theta - factor)
    angle = vifelse(alive, (theta - factor)/2, theta_0)
    s, c = sincos(angle)
    q2 = (c, theta_0, s, theta_0)
    q3 = quat_mul(q2, q1)
    q4 = quat_mul(q3, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q4
  end
end


@makekernel function exact_bend_with_rotation!(i, coords::Coords, e1, e2, theta, a, g, Kn0, w, w_inv, tilde_m, beta_0, L)
  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end
  fringe!(i, coords, a, tilde_m, Kn0, nothing, nothing, e1, e2, 1)
  exact_bend!(i, coords, theta, g, Kn0, tilde_m, beta_0, a, L)
  fringe!(i, coords, a, tilde_m, Kn0, nothing, nothing, e1, e2, -1)
  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end
end