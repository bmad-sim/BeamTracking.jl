"""
Converts the transverse momenta from canonical momenta (px, py) to s derivatives
(x', y').
"""
@makekernel fastgtpsa=true function canonical_to_prime!(i, coords::Coords, g, ax, ay)
  v = coords.v

  rel_p = 1 + v[i,PZI]
  px = v[i,PXI] - ax
  py = v[i,PYI] - ay

  pl2 = rel_p*rel_p - px*px - py*py
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  h = (1 + g*v[i,XI])/pl

  new_px = h*px 
  new_py = h*py

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end


"""
Converts the transverse momenta from s derivatives (x', y') to canonical momenta
(px, py).
"""
@makekernel fastgtpsa=true function prime_to_canonical!(i, coords::Coords, g, ax, ay)
  v = coords.v

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]

  pl2 = h*h + v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  new_px = rel_p*v[i,PXI]/pl + ax
  new_py = rel_p*v[i,PYI]/pl + ay

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end


""" 
Returns the perpendicular component of e_vec divided by the speed of light plus
the cross product of beta and b_vec.
"""
@inline function radiation_field(e_vec, b_vec, beta)
  @inbounds begin @FastGTPSA begin
    e_dot_beta = e_vec[1]*beta[1] + e_vec[2]*beta[2] + e_vec[3]*beta[3]
    e_perp_x = e_vec[1] - e_dot_beta*beta[1]
    e_perp_y = e_vec[2] - e_dot_beta*beta[2]
    e_perp_z = e_vec[3] - e_dot_beta*beta[3]

    beta_cross_b_x = beta[2]*b_vec[3] - beta[3]*b_vec[2]
    beta_cross_b_y = beta[3]*b_vec[1] - beta[1]*b_vec[3]
    beta_cross_b_z = beta[1]*b_vec[2] - beta[2]*b_vec[1]

    field_x = e_perp_x/C_LIGHT + beta_cross_b_x
    field_y = e_perp_y/C_LIGHT + beta_cross_b_y
    field_z = e_perp_z/C_LIGHT + beta_cross_b_z

    return (field_x, field_y, field_z)
  end end
end


"""
Gives radiation damping kick in an electromagnetic field. It is assumed that
the coordinate system has already been rotated such that the curvature
is in the horizontal plane.
"""
@makekernel fastgtpsa=true function deterministic_radiation_field!(i, coords::Coords, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)
  v = coords.v
  canonical_to_prime!(i, coords, g, ax, ay)

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]

  pl2 = h*h + v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  beta = (v[i,PXI]/pl, v[i,PYI]/pl, h/pl)
  field = radiation_field(e_vec, b_vec, beta)
  field_2 = field[1]*field[1] + field[2]*field[2] + field[3]*field[3]

  K = -pl * E_CHARGE/(4*pi*EPS_0) * 2/3 * (q*q)/(mc2*mc2*mc2*mc2) * (E_ref*E_ref*E_ref) * field_2 * L
  new_pz = (v[i,PZI] + rel_p*K)/(1 - rel_p*K)
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])

  prime_to_canonical!(i, coords, g, ax, ay)
end


"""
Gives radiation damping kick in a multipole. It is assumed that
the coordinate system has already been rotated such that the curvature
is in the horizontal plane.
"""
@makekernel fastgtpsa=true function deterministic_radiation_multipole!(i, coords::Coords, q, mc2, E_ref, g, mm, kn, ks, L) 
  v = coords.v

  # Vector potential is (ax, ay, does-not-matter)
  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  zero_0 = zero(kn[1])
  if mm[1] == 0
    b_vec = (bx, by, kn[1])
  else
    b_vec = (bx, by, zero_0)
  end

  e_vec = (zero_0, zero_0, zero_0)  # No electric multipole component

  deterministic_radiation_field!(i, coords, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)
end


"""
Gives radiation diffusion kick in an electromagnetic field. It is assumed that
the coordinate system has already been rotated such that the curvature
is in the horizontal plane.
"""
@makekernel function stochastic_radiation_field!(i, coords::Coords, backend, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)
  v = coords.v

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]
  gamma = rel_p*E_ref/mc2
  px = v[i,PXI] - ax
  py = v[i,PYI] - ay

  pl2 = rel_p*rel_p - px*px - py*py
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  beta = (px/rel_p, py/rel_p, pl/rel_p)
  field = radiation_field(e_vec, b_vec, beta)
  field_2 = field[1]*field[1] + field[2]*field[2] + field[3]*field[3]
  field_1 = sqrt(field_2)
  field_3 = field_2*field_1

  beta_cross_field_x = beta[2]*field[3] - beta[3]*field[2]
  beta_cross_field_y = beta[3]*field[1] - beta[1]*field[3]
  beta_cross_field_z = beta[1]*field[2] - beta[2]*field[1]
  beta_cross_field_2 = beta_cross_field_x*beta_cross_field_x + beta_cross_field_y*beta_cross_field_y + beta_cross_field_z*beta_cross_field_z
  beta_cross_field_1 = sqrt(beta_cross_field_2)

  beta_cross_field_hat_x = vifelse(beta_cross_field_1 > 0, beta_cross_field_x/beta_cross_field_1, 0)
  beta_cross_field_hat_y = vifelse(beta_cross_field_1 > 0, beta_cross_field_y/beta_cross_field_1, 0)

  dt_ds = h*rel_p/pl
  coeff = 55/(24*sqrt(3))/(4*pi*EPS_0)*H_BAR*C_LIGHT*E_CHARGE # H_BAR in eV*s

  mc27 = mc2*mc2*mc2*mc2*mc2*mc2*mc2
  E_ref5 = E_ref*E_ref*E_ref*E_ref*E_ref
  rel_p4 = rel_p*rel_p*rel_p*rel_p
  q2 = q*q

  sigma2 = dt_ds * coeff * q2/mc27 * E_ref5 * rel_p4 * field_3 * L
  sigma2_1 = one(sigma2)
  sigma = sqrt(vifelse(alive, sigma2, sigma2_1)) 

  dpz, theta = gaussian_random(backend, sigma, sqrt(13/55)/gamma)
  s, c = sincos(theta)

  new_px = v[i,PXI] + dpz * (c*beta[1] + s*beta_cross_field_hat_x)
  new_py = v[i,PYI] + dpz * (c*beta[2] + s*beta_cross_field_hat_y)
  new_pz = v[i,PZI] + dpz

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end



"""
Gives radiation diffusion kick in a multipole.
"""
@makekernel fastgtpsa=true function stochastic_radiation!(i, coords::Coords, s, backend, q, mc2, E_ref, g, tilt_ref, mm, kn, ks, L) 
  v = coords.v

  w = rot_quaternion(0, 0, -tilt_ref)
  rotation!(i, coords, w, 0)

  # Vector potential is (ax, ay, does-not-matter)
  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  zero_0 = zero(kn[1])
  if mm[1] == 0
    b_vec = (bx, by, kn[1])
  else
    b_vec = (bx, by, zero_0)
  end

  e_vec = (zero_0, zero_0, zero_0)  # No electric multipole component

  stochastic_radiation_field!(i, coords, backend, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)

  w_inv = inv_rot_quaternion(0, 0, -tilt_ref)
  rotation!(i, coords, w_inv, 0)
end