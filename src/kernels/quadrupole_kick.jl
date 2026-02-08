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
@makekernel fastgtpsa=true function mkm_quadrupole!(i, coords::Coords, q, mc2, radiation_damping, beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  rel_p = 1 + coords.v[i,PZI]
  px = coords.v[i,PXI]
  py = coords.v[i,PYI]
  P_s2 = rel_p*rel_p - px*px - py*py
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  #println(kn)

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, 0, tilde_m, mm, kn, ks, L / 2)
  end

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  multipole_kick!(i, coords, mm, knl, ksl, 2)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  rotation!( i, coords, w, 0)
  quadrupole_matrix!(           i, coords, k1, L)
  rotation!( i, coords, w_inv, 0)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  multipole_kick!(i, coords, mm, knl, ksl, 2)

  if radiation_damping
    deterministic_radiation!(   i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, 0, tilde_m, mm, kn, ks, L / 2)
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

  focus = k1 >= 0  # horizontally focusing if positive

  rel_p = 1 + v[i,PZI]
  xp = v[i,PXI] / rel_p  # x'
  yp = v[i,PYI] / rel_p  # y'
  sqrtks = sqrt(abs(k1 / rel_p)) * s  # |κ|s

  cosine = cos(sqrtks)
  coshine = cosh(sqrtks)
  sinecu = sincu(sqrtks)
  shinecu = sinhcu(sqrtks)
  cx = vifelse(focus, cosine, coshine)
  cy = vifelse(focus, coshine, cosine)
  sx = vifelse(focus, sinecu, shinecu)
  sy = vifelse(focus, shinecu, sinecu)

  new_px = v[i,PXI] * cx - k1 * s * v[i,XI] * sx
  new_py = v[i,PYI] * cy + k1 * s * v[i,YI] * sy
  new_z = v[i,ZI]  - (s / 4) * (  xp*xp * (1 + sx * cx)
                                    + yp*yp * (1 + sy * cy)
                                    + k1 / (1 + v[i,PZI])
                                        * ( v[i,XI]*v[i,XI] * (1 - sx * cx)
                                          - v[i,YI]*v[i,YI] * (1 - sy * cy) )
                                  ) + sign(k1) * ( v[i,XI] * xp * (sqrtks * sx)* 
                                  (sqrtks * sx) - v[i,YI] * yp * (sqrtks * sy)*
                                  (sqrtks * sy) ) / 2
  new_x = v[i,XI] * cx + xp * s * sx
  new_y = v[i,YI] * cy + yp * s * sy
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
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