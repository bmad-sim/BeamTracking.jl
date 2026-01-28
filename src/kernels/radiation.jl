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


@makekernel fastgtpsa=true function deterministic_radiation!(i, coords::Coords, q, mc2, E0, g, mm, kn, ks, L)
  v = coords.v

  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

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

  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  bz_0 = zero(kn[1])
  bz = mm[1] == 0 ? kn[1] : bz_0

  betax = v[i,PXI] / pl
  betay = v[i,PYI] / pl
  betaz = h / pl

  dot = bx*betax + by*betay + bz*betaz

  b_perp_x = bx - dot*betax
  b_perp_y = by - dot*betay
  b_perp_z = bz - dot*betaz

  b_perp_2 = b_perp_x*b_perp_x + b_perp_y*b_perp_y + b_perp_z*b_perp_z

  coeff = E_CHARGE/(4*pi*EPS_0)

  K = -pl * coeff * 2/3 * (q*q)/(mc2*mc2*mc2*mc2) * (E0*E0*E0) * b_perp_2 * L

  new_pz = (v[i,PZI] + rel_p*K)/(1 - rel_p*K)
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])

  prime_to_canonical!(i, coords, g, ax, ay)
end


@makekernel fastgtpsa=true function stochastic_radiation!(i, coords::Coords, q, mc2, E0, g, tilt_ref, mm, kn, ks, L)
  v = coords.v

  w = rot_quaternion(0, 0, tilt_ref)
  rotation!(i, coords, w, 0)

  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  h = 1 + g*v[i,XI]
  rel_p = 1 + v[i,PZI]
  gamma = rel_p * E0 / mc2
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

  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  bz_0 = zero(kn[1])
  bz = mm[1] == 0 ? kn[1] : bz_0

  betax = px / rel_p
  betay = py / rel_p
  betaz = pl / rel_p

  dot = bx*betax + by*betay + bz*betaz

  b_perp_x = bx - dot*betax
  b_perp_y = by - dot*betay
  b_perp_z = bz - dot*betaz

  b_perp_2 = b_perp_x*b_perp_x + b_perp_y*b_perp_y + b_perp_z*b_perp_z
  b_perp = sqrt(b_perp_2)

  dt_ds = h * rel_p / pl

  coeff = 55/(24*sqrt(3))/(4*pi*EPS_0)*H_BAR*C_LIGHT

  mc27 = mc2*mc2*mc2*mc2*mc2*mc2*mc2
  E05 = E0*E0*E0*E0*E0
  rel_p4 = rel_p*rel_p*rel_p*rel_p
  b_perp_3 = b_perp_2*b_perp
  q2 = q*q

  sigma2 = dt_ds * coeff * q2/mc27 * E05 * rel_p4 * b_perp_3 * L
  sigma2_1 = one(sigma2)
  sigma = sqrt(vifelse(alive, sigma2, sigma2_1)) 

  dpz   = gaussian_random(sigma)
  theta = gaussian_random(sqrt(13/55)/gamma)
  s, c  = sincos(theta)

  b_perp_hat_x = vifelse(b_perp > 0, b_perp_x / b_perp, 0)
  b_perp_hat_y = vifelse(b_perp > 0, b_perp_y / b_perp, 0)

  new_px = v[i,PXI] + dpz * (c*betax + s*b_perp_hat_x)
  new_py = v[i,PYI] + dpz * (c*betay + s*b_perp_hat_y)
  new_pz = v[i,PZI] + dpz

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])

  w_inv = inv_rot_quaternion(0, 0, tilt_ref)
  rotation!(i, coords, w_inv, 0)
end