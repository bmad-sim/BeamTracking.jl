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
@makekernel fastgtpsa=true function bkb_multipole!(i, coords::Coords, q, mc2, radiation_damping, tilde_m, beta_0, a, g, w, w_inv, k0, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  rotation!( i, coords, w, 0)

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, g, tilde_m, mm, kn, ks, L / 2)
  end

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  end

  multipole_kick!(i, coords, mm, knl, ksl, 1)
  exact_bend!(    i, coords, g*L, g, k0, tilde_m, beta_0, L)
  multipole_kick!(i, coords, mm, knl, ksl, 1)

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, g, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, g, tilde_m, mm, kn, ks, L / 2)
  end

  rotation!( i, coords, w_inv, 0)
end 

"""
    exact_bend!(i, coords::Coords, e1, e2, theta, g, Kn0, w, w_inv, tilde_m, beta_0, L)

Tracks a particle through a sector bend via exact tracking. If edge angles are 
provided, a linear hard-edge fringe map is applied at both ends.

#Arguments
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'Kn0'      -- normalized dipole field
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'L'        -- length
"""
@makekernel fastgtpsa=true function exact_bend!(i, coords::Coords, theta, g, Kn0, tilde_m, beta_0, L)
  v = coords.v
  rel_p = 1 + v[i,PZI]

  pt2 = rel_p*rel_p - v[i,PYI]*v[i,PYI]
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
  cplus = cos(phi1) 
  splus = sin(phi1)
  sinc_theta = sincu(theta)
  sinc_theta_2 = sincu(theta/2)
  cosc_theta = sinc_theta_2*sinc_theta_2/2
  sgn = sign(L)
  alpha_helper = h*L*sinc_theta
  alpha = 2*h*splus*L*sinc_theta - gp*alpha_helper*alpha_helper

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
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
end


@makekernel fastgtpsa=true function linear_bend_fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, e, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  f = Kn0*tan(e)

  if !isnothing(coords.q)
    b_vec = (-v[i,YI]*f - sign*v[i,XI]*Ksol/2, -v[i,XI]*f - sign*v[i,YI]*Ksol/2, sign*v[i,YI]*Kn0)
    ax = -v[i,YI]*Ksol/2
    ay =  v[i,XI]*Ksol/2
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
  end

  new_px = v[i,PXI] + f*v[i,XI]
  new_py = v[i,PYI] - f*v[i,YI]
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])

  if !isnothing(coords.q)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
  end
end


@makekernel fastgtpsa=true function exact_bend_with_rotation!(i, coords::Coords, e1, e2, theta, a, g, Kn0, w, w_inv, tilde_m, beta_0, L)
  rotation!(i, coords, w, 0)
  linear_bend_fringe!(i, coords, a, tilde_m, 0, Kn0, e1, 1)
  exact_bend!(i, coords, theta, g, Kn0, tilde_m, beta_0, L)
  linear_bend_fringe!(i, coords, a, tilde_m, 0, Kn0, e2, -1)
  rotation!(i, coords, w_inv, 0)
end


# This is separate because the spin can be transported exactly here
@makekernel fastgtpsa=true function exact_curved_drift!(i, coords::Coords, e1, e2, theta, g, w, w_inv, a, tilde_m, beta_0, L) 
  exact_bend_with_rotation!(i, coords, 0, 0, theta, a, g, 0, w, w_inv, tilde_m, beta_0, L)
  if !isnothing(coords.q)
    rotation!(i, coords, w, 0)
    rotate_spin!(i, coords, a, g, tilde_m, SA[0], SA[0], SA[0], L)
    rotation!(i, coords, w_inv, 0)
  end
end
