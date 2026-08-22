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
  new_tau = v[i,ZI] + L/beta_0 - (x*sinh_arg + ptau_plus_inv_beta_0*sinh_arg_over_kE + px*cosh_arg_minus_1_over_kE) # make stable

  v[i,XI]  = vifelse(alive, new_x,   v[i,XI])
  v[i,PXI] = vifelse(alive, new_px,  v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,   v[i,YI])
  v[i,ZI]  = vifelse(alive, new_tau, v[i,ZI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    dpx = px*kE*cosh_arg_minus_1_over_kE + (x*kE + ptau_plus_inv_beta_0)*sinh_arg
    p_perp2 = py2 + ps2
    p_perp = sqrt(p_perp2) # protect
    factor = tilde_m2 + py2 + ps2
    mgamma_i = sqrt(factor + px2) # protect
    mgamma_f = sqrt(factor + v[i,PXI]*v[i,PXI]) # protect
    mdgamma = (px + v[i,PXI])/(mgamma_i + mgamma_f)*dpx

    I1 = log((v[i,PXI] + mgamma_f)/(px + mgamma_i))/tilde_m # protect
    I2 = 2*atan(p_perp*(dpx + mdgamma)/(p_perp2 + (px + mgamma_i + tilde_m)*(v[i,PXI] + mgamma_f + tilde_m)))/p_perp # protect
    coeff = a*I1 + I2

    o1 =  zero(coeff)
    o2 =  coeff*ps
    o3 = -coeff*py

    q1 = expq((o1, o2, o3), alive)
    q2 = quat_mul(q1, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q2
  end
end