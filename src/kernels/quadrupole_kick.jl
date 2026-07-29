"""
mkm_quadrupole!()

This integrator uses Matrix-Kick-Matrix to implement a quadrupole
integrator accurate though second-order in the integration step-size. The vectors
kn and ks contain the normal and skew multipole strengths.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
ks: vector of skew multipole strengths scaled by Bρ0
L: element length
"""
@makekernel fastgtpsa=true function mkm_quadrupole!(i, coords::Coords, s, radiation_params, beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L)
  other_multipoles = (length(mm) > 1)
  knl = kn .* L ./ 2
  ksl = ks .* L ./ 2
  
  rel_p = 1 + coords.v[i,PZI]
  px = coords.v[i,PXI]
  py = coords.v[i,PYI]
  P_s2 = rel_p*rel_p - px*px - py*py
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])

  if !isnothing(coords.q) && other_multipoles
    rotate_spin!(i, coords, a, 0, tilde_m, mm, knl, ksl, 2)
  end

  if !isnothing(radiation_params)
    q, mc2, E_ref = radiation_params
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
  end

  if other_multipoles
    multipole_kick!(i, coords, mm, knl, ksl, 2)
  end

  quadrupole_kick!(i, coords, beta_0, gamsqr_0, tilde_m, L / 2)

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  if !isnothing(coords.q)
    quadrupole_magnus6!(i, coords, k1, tilde_m, a, L)
  else
    quadrupole_matrix!(i, coords, k1, L)
  end

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end

  quadrupole_kick!(i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  
  if other_multipoles
    multipole_kick!(i, coords, mm, knl, ksl, 2)
  end

  if !isnothing(radiation_params)
    deterministic_radiation_multipole!(i, coords, q, mc2, E_ref, 0, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q) && other_multipoles
    rotate_spin!(i, coords, a, 0, tilde_m, mm, knl, ksl, 2)
  end
end


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k1:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic p_over_q_ref.
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_matrix!(i, coords::Coords, k1, s)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  focus = (k1 >= 0)  # horizontally focusing if positive

  rel_p = 1 + v[i,PZI]
  xp = v[i,PXI] / rel_p  # x'
  yp = v[i,PYI] / rel_p  # y'
  arg = k1*s*s/rel_p

  sinecu,  cosine  = sincos_quaternion( arg)
  shinecu, coshine = sincos_quaternion(-arg)
  cx = vifelse(focus, cosine,  coshine)
  cy = vifelse(focus, coshine, cosine)
  sx = vifelse(focus, sinecu,  shinecu)
  sy = vifelse(focus, shinecu, sinecu)

  new_px = v[i,PXI] * cx - k1 * v[i,XI] * s * sx
  new_py = v[i,PYI] * cy + k1 * v[i,YI] * s * sy
  new_z  = v[i,ZI]  - (s / 4) * (  xp * xp * (1 + sx * cx)
                                    + yp * yp * (1 + sy * cy)
                                    + k1 / rel_p
                                        * ( v[i,XI]*v[i,XI] * (1 - sx * cx)
                                          - v[i,YI]*v[i,YI] * (1 - sy * cy) )
                                  ) + arg * (v[i,XI] * xp * sx * sx
                                  - v[i,YI] * yp * sy * sy) / 2
  new_x  = v[i,XI] * cx + xp * s * sx
  new_y  = v[i,YI] * cy + yp * s * sy
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])
  v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
  v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
end 


"""
quadrupole_kick!()

Track "remaining part" of quadrupole —— a position kick.

### Note re implementation:
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_kick!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, s)
  v = coords.v

  P      = 1 + v[i,PZI]             # [scaled] total momentum, P/P0 = 1 + δ
  PtSqr  = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]  # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps2    = P*P - PtSqr        
  good_momenta = (Ps2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  Ps2_1 = one(Ps2)
  Ps = sqrt(vifelse(good_momenta, Ps2, Ps2_1)) # longitudinal momentum,   Ps   = √[(1 + δ)^2 - P⟂^2]
  alive = (coords.state[i] == STATE_ALIVE)

  new_x = v[i,XI] + s * v[i,PXI] * PtSqr / (P * Ps * (P + Ps))
  new_y = v[i,YI] + s * v[i,PYI] * PtSqr / (P * Ps * (P + Ps))
  new_z = v[i,ZI] - s * (P * (PtSqr - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(P*P + tilde_m*tilde_m) * Ps
                                    * (beta_0 * sqrt(P*P + tilde_m*tilde_m) + Ps)
                                  )
                            - PtSqr / (2 * P*P))
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end 


# There seems to be a bug with FastGTPSA temps not being released
@makekernel fastgtpsa=false function quadrupole_magnus6!(i, coords::Coords, k1, tilde_m, a, L)
  v = coords.v
  q = coords.q

  rel_p = 1 + v[i,PZI]
  rel_p2 = rel_p*rel_p

  eta = sqrt(15)/10*L

  quadrupole_matrix!(i, coords, k1, L/2 - eta)

  x_1  = v[i,XI]
  px_1 = v[i,PXI]
  y_1  = v[i,YI]
  py_1 = v[i,PYI]

  ps_1_2 = rel_p2 - px_1*px_1 - py_1*py_1
  good_momenta = (ps_1_2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_1_1 = one(ps_1_2)
  ps_1 = sqrt(vifelse(good_momenta, ps_1_2, ps_1_1))

  quadrupole_matrix!(i, coords, k1, eta)

  x_2  = v[i,XI]
  px_2 = v[i,PXI]
  y_2  = v[i,YI]
  py_2 = v[i,PYI]

  ps_2_2 = rel_p2 - px_2*px_2 - py_2*py_2
  good_momenta = (ps_2_2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_2_1 = one(ps_2_2)
  ps_2 = sqrt(vifelse(good_momenta, ps_2_2, ps_2_1))

  quadrupole_matrix!(i, coords, k1, eta)

  x_3  = v[i,XI]
  px_3 = v[i,PXI]
  y_3  = v[i,YI]
  py_3 = v[i,PYI]

  ps_3_2 = rel_p2 - px_3*px_3 - py_3*py_3
  good_momenta = (ps_3_2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_3_1 = one(ps_3_2)
  ps_3 = sqrt(vifelse(good_momenta, ps_3_2, ps_3_1))

  quadrupole_matrix!(i, coords, k1, L/2 - eta)

  beta_gamma = rel_p/tilde_m
  beta_gamma2 = beta_gamma*beta_gamma
  gamma_minus_1 = beta_gamma2/(1 + sqrt(1 + beta_gamma2))
  gamma = gamma_minus_1 + 1
  chi = 1 + a*gamma
  k1_chi = -k1*chi
  k1_g = k1*a*gamma_minus_1/rel_p2

  coeff1_1 = k1_chi/ps_1
  coeff1_2 = k1_chi/ps_2
  coeff1_3 = k1_chi/ps_3

  coeff2_1 = k1_g*(x_1*py_1 + y_1*px_1)/ps_1
  coeff2_2 = k1_g*(x_2*py_2 + y_2*px_2)/ps_2
  coeff2_3 = k1_g*(x_3*py_3 + y_3*px_3)/ps_3

  a1 = (coeff1_1*y_1 + coeff2_1*px_1, coeff1_1*x_1 + coeff2_1*py_1, coeff2_1*ps_1)
  a2 = (coeff1_2*y_2 + coeff2_2*px_2, coeff1_2*x_2 + coeff2_2*py_2, coeff2_2*ps_2)
  a3 = (coeff1_3*y_3 + coeff2_3*px_3, coeff1_3*x_3 + coeff2_3*py_3, coeff2_3*ps_3)

  alpha1 = L .* a2
  alpha2 = (eta*10/3) .* (a3 .- a1)
  alpha3 = (10*L/3) .* (a3 .- (2 .* a2) .+ a1)

  c1 = cross(alpha1, alpha2)
  c2 = cross(alpha1, (2 .* alpha3) .+ c1) ./ (-60)

  omega = alpha1 .+ (alpha3 ./ 12) .+ (cross((-20 .* alpha1) .- alpha3 .+ c1, alpha2 .+ c2) ./ 240)

  q1 = expq(omega, alive)
  q2 = quat_mul(q1, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
  q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q2
end


function cross(v1, v2)
  @inbounds begin @FastGTPSA begin
    a1, b1, c1 = v1
    a2, b2, c2 = v2
    o1 = b1*c2 - b2*c1
    o2 = c1*a2 - c2*a1
    o3 = a1*b2 - a2*b1
  end end
  return (o1, o2, o3)
end