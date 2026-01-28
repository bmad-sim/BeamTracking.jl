
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
@makekernel fastgtpsa=true function sks_multipole!(i, coords::Coords, q, mc2, radiation_damping, beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  exact_solenoid!(  i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)

  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    multipole_kick!(i, coords, mm, knl, ksl, -1)
    rotate_spin!(                 i, coords, a, 0, tilde_m, mm, kn, ks, L)
    multipole_kick!(i, coords, mm, knl, ksl, -1)
  end

  if radiation_damping
    deterministic_radiation!(     i, coords, q, mc2, E0, 0, mm, kn, ks, L / 2)
  end

  exact_solenoid!(  i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
end 

@makekernel fastgtpsa=true function exact_solenoid!(i, coords::Coords, ks, beta_0, gamsqr_0, tilde_m, L)
  v = coords.v

  # Recurring variables
  rel_p = 1 + v[i,PZI]
  px_k = v[i,PXI] + v[i,YI] * ks / 2
  py_k = v[i,PYI] - v[i,XI] * ks / 2
  pt2 = px_k*px_k + py_k*py_k
  pr2 = rel_p*rel_p - pt2
  good_momenta = (pr2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pr2_1 = one(pr2)
  pr = sqrt(vifelse(good_momenta, pr2, pr2_1))

  s = sin(ks * L / pr)
  cp = 1 + cos(ks * L / pr)
  cm = 2 - cp
  # Temporaries
  x_0 = v[i,XI]
  px_0 = v[i,PXI]
  y_0 = v[i,YI]

  new_z = v[i,ZI] - rel_p * L * (pt2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) /
                ( beta_0 * sqrt(rel_p*rel_p + tilde_m*tilde_m) * pr * (beta_0 * 
                sqrt(rel_p*rel_p + tilde_m*tilde_m) + pr) )
  new_x = cp * x_0 / 2 + s * (px_0 / ks + y_0 / 2) + cm * v[i,PYI] / ks
  new_px = s * (v[i,PYI] / 2 - ks * x_0 / 4) + cp * px_0 / 2 - ks * cm * y_0 / 4
  new_y = s * (v[i,PYI] / ks - x_0 / 2) + cp * y_0 / 2 - cm * px_0 / ks
  new_py = ks * cm * x_0 / 4 - s * (px_0 / 2 + ks * y_0 / 4) + cp * v[i,PYI] / 2
  # Update
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end