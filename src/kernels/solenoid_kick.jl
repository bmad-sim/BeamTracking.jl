
"""
sks_multipole!()

This integrator uses Solenoid-Kick-Solenoid to track a beam through
a straight, finite-length multipole magnet with a solenoid field. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the dipole component.

Arguments
—————————
Ksol: solenoid strength
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@makekernel fastgtpsa=true function sks_multipole!(i, coords::Coords, s, radiation_params, beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, ks, L)
  knL = kn .* (L / 2)
  ksL = ks .* (L / 2)

  if !isnothing(radiation_params)
    q, mc2, E_ref = radiation_params
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, kn, ks, 0, L / 2)
  end

  multipole_kick!(i, coords, mm, knL, ksL, 0)
  exact_solenoid!(i, coords, Ksol, beta_0, gamsqr_0, tilde_m, a, L)
  multipole_kick!(i, coords, mm, knL, ksL, 0)

  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, kn, ks, 0, L / 2)
  end

  if !isnothing(radiation_params)
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
  end
end 


@makekernel fastgtpsa=true function exact_solenoid!(i, coords::Coords, ks, beta_0, gamsqr_0, tilde_m, a, L)
  v = coords.v

  # Recurring variables
  rel_p = 1 + v[i,PZI]
  rel_p2 = rel_p*rel_p
  px_k = v[i,PXI] + v[i,YI] * ks / 2
  py_k = v[i,PYI] - v[i,XI] * ks / 2
  pt2 = px_k*px_k + py_k*py_k
  pr2 = rel_p2 - pt2
  good_momenta = (pr2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pr2_1 = one(pr2)
  pr = sqrt(vifelse(good_momenta, pr2, pr2_1))

  arg = ks*L/pr
  s_over_ks = L/pr*sincu(arg)
  s = ks*s_over_ks
  cm_over_ks = arg*L/pr*one_cos_norm(arg)
  cm_times_ks = ks*ks*cm_over_ks
  cp = 2 - ks*cm_over_ks
  # Temporaries
  x_0 = v[i,XI]
  px_0 = v[i,PXI]
  y_0 = v[i,YI]

  new_z = v[i,ZI] - rel_p * L * (pt2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) /
                ( beta_0 * sqrt(rel_p2 + tilde_m*tilde_m) * pr * (beta_0 * 
                sqrt(rel_p2 + tilde_m*tilde_m) + pr) )
  new_x = cp * x_0 / 2 + s_over_ks * (px_0 + y_0 * ks / 2) + cm_over_ks * v[i,PYI]
  new_px = s * (v[i,PYI] / 2 - x_0 * ks / 4) + cp * px_0 / 2 - cm_times_ks * y_0 / 4
  new_y = s_over_ks * (v[i,PYI] - x_0 * ks / 2) + cp * y_0 / 2 - cm_over_ks * px_0
  new_py = cm_times_ks * x_0 / 4 - s * (px_0 / 2 + y_0 * ks / 4) + cp * v[i,PYI] / 2
  # Update
  v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])
  v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    beta_gamma = rel_p/tilde_m
    beta_gamma2 = beta_gamma*beta_gamma
    gamma_minus_1 = beta_gamma2/(1 + sqrt(1 + beta_gamma2))
    coeff = a*ks*L*gamma_minus_1/rel_p2

    o1 = coeff*px_k
    o2 = coeff*py_k
    o3 = -a*arg*(1 + pt2*gamma_minus_1/rel_p2)
    q1 = expq((o1, o2, o3), alive)

    arg_0 = zero(arg)
    angle = vifelse(alive, -arg/2, arg_0)
    s, c = sincos(angle)
    q2 = (c, arg_0, arg_0, s)
    q3 = quat_mul(q2, q1)
    q4 = quat_mul(q3, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q4
  end
end