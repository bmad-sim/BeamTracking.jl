"""
This function computes the integrated spin-precession vector using the multipole 
coefficients kn and ks indexed by mm, i.e., knl[i] is the normal 
coefficient of order mm[i].
"""
function omega_multipole(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v

    if mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = zero(v[i,XI])
      ay = ax
    end

    bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
    bz_0 = zero(kn[1])
    bz = mm[1] == 0 ? kn[1] : bz_0
    b_vec = (bx, by, bz)
    e_vec = (bz_0, bz_0, bz_0)

    omega = omega_field(i, coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
  end end

  return omega
end


"""
This function computes the integrated spin-precession vector using the fields
and returns it together with the endpoint longitudinal-momentum validity flag.
"""
@inline function omega_field_core(x, px, py, pz, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
  pkin_x = px - ax
  pkin_y = py - ay
  rel_p = 1 + pz
  beta_gamma = rel_p / tilde_m
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma / gamma

  pl2 = rel_p*rel_p - pkin_x*pkin_x - pkin_y*pkin_y
  good_momenta = (pl2 > zero(pl2))
  pl = sqrt(vifelse(good_momenta, pl2, one(pl2)))

  coeff = -(1 + g*x)/pl
  coeff1 = coeff * (1 + a*gamma)
  coeff2 = coeff * (1 + a)
  coeff3 = -coeff * beta / C_LIGHT * gamma * (a + 1/(1+gamma))/rel_p

  betax = pkin_x / rel_p
  betay = pkin_y / rel_p
  betaz = pl / rel_p

  dot = b_vec[1]*betax + b_vec[2]*betay + b_vec[3]*betaz

  b_para_x = dot * betax
  b_para_y = dot * betay
  b_para_z = dot * betaz

  b_perp_x = (b_vec[1] - b_para_x) * coeff1
  b_perp_y = (b_vec[2] - b_para_y) * coeff1
  b_perp_z = (b_vec[3] - b_para_z) * coeff1

  b_para_x = b_para_x * coeff2
  b_para_y = b_para_y * coeff2
  b_para_z = b_para_z * coeff2

  e_part_x = (pkin_y*e_vec[3] - pl*e_vec[2]) * coeff3
  e_part_y = (pl*e_vec[1] - pkin_x*e_vec[3]) * coeff3
  e_part_z = (pkin_x*e_vec[2] - pkin_y*e_vec[1]) * coeff3

  ox = (b_perp_x + b_para_x + e_part_x) * L
  oy = (b_perp_y + b_para_y + e_part_y + g) * L
  oz = (b_perp_z + b_para_z + e_part_z) * L

  return (ox, oy, oz), good_momenta
end

function omega_field(i, coords::Coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    omega, good_momenta = omega_field_core(
      v[i,XI], v[i,PXI], v[i,PYI], v[i,PZI],
      a, g, tilde_m, ax, ay, e_vec, b_vec, L
    )
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  end end
  return omega
end

"""
This function rotates particle i's quaternion according to the multipoles present.
"""
@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_multipole(i, coords, a, g, tilde_m, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end


@makekernel fastgtpsa=true function integrate_with_spin_thin!(i, coords::Coords, ker, params, a, g, tilde_m, mm, knl, ksl)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  ker(i, coords, params...)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
end


@makekernel fastgtpsa=true function rotate_spin_field!(i, coords::Coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_field(i, coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end