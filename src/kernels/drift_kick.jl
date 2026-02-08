"""
dkd_multipole!()

This integrator uses Drift-Kick-Drift to track a beam through
a straight, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the dipole component. (For example, kn[3] denotes
the normal sextupole strength in Tesla/m^2.) The argument ns
denotes the number of slices.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
ks: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@makekernel fastgtpsa=true function dkd_multipole!(i, coords::Coords, q, mc2, radiation_damping, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  E0 = mc2/tilde_m/beta_0 # could probably exclude beta_0 because ultrarelativistic radiation

  exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)

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

  exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
end

"""
    exact_drift!(i, coords, β_0, γsqr_0, tilde_m, L)

Return the result of exact tracking a particle through a drift
of length `L`, assuming `β_0`, `γsqr_0`, and `tilde_m` respectively
denote the reference velocity normalized to the speed of light,
the corresponding value of the squared Lorentz factor, and the
particle rest energy normalized to the reference value of ``pc``.

NB: In the computation of ``z_final``, we use the fact that
  - ``1/√a - 1/√b == (coords - a)/(√a √b (√a + √b))``
to avoid the potential for severe cancellation when
``a`` and ``coords`` both have the form ``1 + ε`` for different small
values of ``ε``.

## Arguments
- `β_0`:     reference velocity normalized to the speed of light, ``v_0 / c``
- `γsqr_0`:  corresponding value of the squared Lorentz factor
- `tilde_m`: particle rest energy normalized to the reference value of ``pc``
- `L`:       element length, in meters
"""
@makekernel fastgtpsa=true function exact_drift!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, L)
  v = coords.v

  rel_p = 1 + v[i,PZI]
  P_t2 = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
  P_s2 = rel_p*rel_p - P_t2
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  P_s2_1 = one(P_s2)
  P_s = sqrt(vifelse(good_momenta, P_s2, P_s2_1))

  new_x = v[i,XI] + v[i,PXI] * L / P_s
  new_y = v[i,YI] + v[i,PYI] * L / P_s
  # high-precision computation of z_final:
  #   vf.z = vi.z - (1 + δ) * L * (1 / Ps - 1 / (β0 * sqrt((1 + δ)^2 + tilde_m^2)))
  new_z = v[i,ZI] - ( (1 + v[i,PZI]) * L
                * (P_t2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                / ( beta_0 * sqrt(rel_p*rel_p + tilde_m*tilde_m) * P_s
                    * (beta_0 * sqrt(rel_p*rel_p + tilde_m*tilde_m) + P_s)
                  )
              )
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end

function drift_params(species::Species, p_over_q_ref)
  beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0*beta_gamma_0
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end
