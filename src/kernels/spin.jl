#---------------------------------------------------------------------------------------------------
"""
This function computes the integrated spin-precession vector using the magnetic multipole 
coefficients kn and ks indexed by mm.
"""
function omega_multipole(i, coords::Coords, a, g, tilde_m, mm, kn, ks, excluding, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v

    # Vector potential is (ax, ay, does-not-matter)
    if mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = zero(v[i,XI])
      ay = ax
    end

    bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], excluding)
    zero_0 = zero(kn[1])
    if mm[1] == 0 && excluding != 0
      b_vec = (bx, by, kn[1])
    else
      b_vec = (bx, by, zero_0)
    end

    e_vec = (zero_0, zero_0, zero_0)   # No electric multipole component

    omega = omega_field(i, coords, a, g, tilde_m, ax, ay, e_vec, b_vec, Val{false}(), L)
  end end

  return omega
end

"""
This function computes the integrated spin-precession vector using the integrated magnetic multipole 
coefficients KnL and KsL indexed by mm.

Any solenoid component is ignored.
"""
function omega_multipole(i, coords::Coords, a, g, tilde_m, mm, KnL, KsL, excluding)
  @FastGTPSA begin @inbounds begin
    v = coords.v

    bx, by = normalized_field(mm, KnL, KsL, v[i,XI], v[i,YI], excluding)
    zero_0 = zero(KnL[1])
    b_vec = (bx, by, zero_0)
    e_vec = (zero_0, zero_0, zero_0)   # No electric multipole component

    omega = omega_field(i, coords, a, g, tilde_m, zero(v[i,XI]), zero(v[i,XI]), e_vec, b_vec, Val{false}(), 1)
  end end

  return omega
end

#---------------------------------------------------------------------------------------------------

"""
    omega_field(i, coords::Coords, a, g, tilde_m, ax, ay, e_vec, b_vec, include_curvature, L) -> omega_vec

This function computes the integrated spin-precession vector using the fields.

## Input:

- `a`                  Anomalous magnetic moment.
- `g`                  Reference bend strength 1/radius when in a bend element.
- `tilde_m             Normalized mass `mass / P0c`
- `ax`, `ay`           Transverse vector potential components.
- `e_vec`              Normalized electric field `e_field * q / P_ref`
- `b_vec`              Normalized magnetic field `b_field * q / P_ref`
- `include_curvature`  If true, include the curvature term in the spin precession.

## Output:
- `omega_vec  3D Rotation tuple. 
"""
function omega_field(i, coords::Coords, a, g, tilde_m, ax, ay, e_vec, b_vec, ::Val{include_curvature}, L) where {include_curvature}
  @FastGTPSA begin @inbounds begin
    v = coords.v
    px = v[i,PXI] - ax
    py = v[i,PYI] - ay
    rel_p = 1 + v[i,PZI]
    beta_gamma = rel_p / tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    beta = beta_gamma / gamma

    pl2 = rel_p*rel_p - px*px - py*py
    pl2_0 = zero(pl2)
    good_momenta = (pl2 > pl2_0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    pl2_1 = one(pl2)
    pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

    coeff = -(1 + g*v[i,XI])/pl
    coeff1 = coeff * (1 + a*gamma)
    coeff2 = coeff * (1 + a)
    coeff3 = -coeff * beta * gamma * (a + 1/(1+gamma))/C_LIGHT

    betax = px / rel_p
    betay = py / rel_p
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

    e_part_x = (betay*e_vec[3] - betaz*e_vec[2]) * coeff3
    e_part_y = (betaz*e_vec[1] - betax*e_vec[3]) * coeff3
    e_part_z = (betax*e_vec[2] - betay*e_vec[1]) * coeff3

    ox = (b_perp_x + b_para_x + e_part_x) * L        
    oy = (b_perp_y + b_para_y + e_part_y + vifelse(include_curvature, g, zero(g))) * L
    oz = (b_perp_z + b_para_z + e_part_z) * L

    omega = (ox, oy, oz)
  end end
  return omega
end

#---------------------------------------------------------------------------------------------------
"""
This function rotates particle i's quaternion according to the multipoles present.
"""
@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, KnL, KsL, excluding)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_multipole(i, coords, a, g, tilde_m, mm, KnL, KsL, excluding), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end

@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, Kn, Ks, excluding, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_multipole(i, coords, a, g, tilde_m, mm, Kn, Ks, excluding, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end

#---------------------------------------------------------------------------------------------------


@makekernel fastgtpsa=true function rotate_spin_field!(i, coords::Coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_field(i, coords, a, g, tilde_m, ax, ay, e_vec, b_vec, Val{true}(), L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end