@makekernel fastgtpsa=true function exact_elsep!(i, coords::Coords, kE, beta_0, tilde_m, a, L)
  v = coords.v
  x = v[i,XI]
  px = v[i,PXI]
  px2 = px*px
  py = v[i,PYI]
  py2 = py*py
  ptau_plus_inv_beta_0 = v[i,PZI] + 1/beta_0
  rel_E = ptau_plus_inv_beta_0 + kE*v[i,XI]
  tilde_m2 = tilde_m*tilde_m
  
  ps2 = rel_E*rel_E - tilde_m2 - px2 - py2
  good_momenta = (ps2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps2_1 = one(ps2)
  ps = sqrt(vifelse(good_momenta, ps2, ps2_1))
  
  arg = kE*L/ps
  sinhcu_arg_half = sinhcu(arg/2)
  sinhcu_arg_half_2 = sinhcu_arg_half*sinhcu_arg_half
  cosh_arg_minus_1_over_kE = kE*L*L/(2*ps2)*sinhcu_arg_half_2
  cosh_arg = 1 + kE*cosh_arg_minus_1_over_kE
  sinh_arg_over_kE = L/ps*sinhcu_arg_half*sqrt(1 + sinhcu_arg_half_2*arg*arg/4)
  sinh_arg = kE*sinh_arg_over_kE

  new_x   = x*cosh_arg + ptau_plus_inv_beta_0*cosh_arg_minus_1_over_kE + px*sinh_arg_over_kE
  new_px  = px*cosh_arg + (x*kE + ptau_plus_inv_beta_0)*sinh_arg
  new_y   = v[i,YI] + L*py/ps
  new_tau = v[i,ZI] + L/beta_0 - (x*sinh_arg + ptau_plus_inv_beta_0*sinh_arg_over_kE + px*cosh_arg_minus_1_over_kE) # probably not optimal

  v[i,XI]  = vifelse(alive, new_x,   v[i,XI])
  v[i,PXI] = vifelse(alive, new_px,  v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,   v[i,YI])
  v[i,ZI]  = vifelse(alive, new_tau, v[i,ZI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    dpx = px*kE*cosh_arg_minus_1_over_kE + (x*kE + ptau_plus_inv_beta_0)*sinh_arg
    p_perp2 = py2 + ps2
    factor = tilde_m2 + p_perp2

    good_momenta = (p_perp2 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    p_perp = sqrt(vifelse(good_momenta, p_perp2, ps2_1))

    mgamma_i2 = factor + px2
    good_momenta = (mgamma_i2 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    mgamma_i = sqrt(vifelse(good_momenta, mgamma_i2, ps2_1))

    mgamma_f2 = factor + v[i,PXI]*v[i,PXI]
    good_momenta = (mgamma_f2 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    mgamma_f = sqrt(vifelse(good_momenta, mgamma_f2, ps2_1))

    mdgamma = (px + v[i,PXI])/(mgamma_i + mgamma_f)*dpx

    arg1 = vifelse(alive, (v[i,PXI] + mgamma_f)/(px + mgamma_i), ps2_1)
    I1 = log(arg1)/tilde_m
    I2 = 2*atan(p_perp*(dpx + mdgamma)/(p_perp2 + (px + mgamma_i + tilde_m)*(v[i,PXI] + mgamma_f + tilde_m)))/p_perp
    coeff = a*I1 + I2

    o1 =  zero(coeff)
    o2 =  coeff*ps
    o3 = -coeff*py

    q1 = expq((o1, o2, o3), alive)
    q2 = quat_mul(q1, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q2
  end
end