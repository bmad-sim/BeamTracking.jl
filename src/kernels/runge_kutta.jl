# Reference factors are constant throughout a particle's RK tracking loop.
@inline function _rk_reference_factors(beta_0, tilde_m)
  inv_beta0 = inv(beta_0)
  inv_c = inv(c_light(typeof(tilde_m)))
  return (; inv_beta0, inv_c)
end

# Factors shared by the field-query time and equations of motion at one RK stage.
@inline function _rk_kinematics(pz, tilde_m)
  rel_p = 1 + pz
  valid_rel_p = (rel_p > 0) & (rel_p < Inf)
  inv_rel_p = inv(vifelse(rel_p > 0, rel_p, one(rel_p)))
  energy_squared = rel_p^2 + tilde_m^2
  energy = sqrt(energy_squared)
  beta = rel_p / energy
  inv_beta = energy * inv_rel_p
  return (; rel_p, valid_rel_p, inv_rel_p, energy_squared, beta, inv_beta)
end

# Element-local time: the reference particle enters at t = 0.
@inline function _rk_time_from_factors(z, s, tilde_m, beta_0, kinematics, reference=_rk_reference_factors(beta_0, tilde_m))
  # Keep rejected particles' query times finite for zero/negative/infinite momentum.
  fallback_inv_beta = sqrt(one(tilde_m) + tilde_m^2)
  time_inv_beta = vifelse(kinematics.valid_rel_p, kinematics.inv_beta, fallback_inv_beta)
  return (s * reference.inv_beta0 - z * time_inv_beta) * reference.inv_c
end

@inline function _rk_particle_time(z, pz, s, tilde_m, beta_0)
  return _rk_time_from_factors(z, s, tilde_m, beta_0, _rk_kinematics(pz, tilde_m))
end

# Mechanical momenta must describe forward motion with nonzero longitudinal momentum.
@inline function _valid_momentum(px, py, pz)
  rel_p = 1 + pz
  inv_p = inv(vifelse(rel_p > 0, rel_p, one(rel_p)))
  return (rel_p > 0) & (rel_p < Inf) & ((px * inv_p)^2 + (py * inv_p)^2 < 1)
end

"""
  _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
        tilde_m, beta_0, gx, gy)

Calculate the derivative vector du/ds for relativistic particle tracking.
Returns an SVector{6} containing [dx/ds, dpx/ds, dy/ds, dpy/ds, dz/ds, dpz/ds].

Uses branchless operations for GPU/SIMD compatibility. For unphysical momenta,
returns zero derivatives (caller should mark particle as lost).

# Arguments
- `x, px, y, py, z, pz`: State vector components
- `s`: Arc length position
- `Ex, Ey, Ez`: Electric field divided by reference rigidity
- `Bx, By, Bz`: Magnetic field divided by reference rigidity
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gx`, `gy`: Horizontal and vertical reference curvature components
"""
@inline function _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, gx, gy, kinematics, reference=_rk_reference_factors(beta_0, tilde_m))
  rel_p = kinematics.rel_p
  inv_rel_p = kinematics.inv_rel_p

  # Transverse velocity components (normalized)
  vt_x = px * inv_rel_p
  vt_y = py * inv_rel_p
  vt2 = vt_x^2 + vt_y^2

  # Check for unphysical momenta (branchless)
  vt2_1 = one(vt2)
  good_momenta = (rel_p > 0) & (rel_p < Inf) & (vt2 < vt2_1)
  vt2_safe = vifelse(good_momenta, vt2, zero(vt2))

  # Particle beta and energy were computed before the field query.
  beta = kinematics.beta
  
  # Unit direction and path-length factor.
  uz = sqrt(1 - vt2_safe)
  dh_bend = x * gx + y * gy
  inv_uz = inv(uz)
  path_factor = (1 + dh_bend) * inv_uz
  inv_beta = kinematics.inv_beta
  electric_factor = inv_beta * reference.inv_c

  dx_ds = vt_x * path_factor
  dy_ds = vt_y * path_factor
  dpx_ds = (Ex * electric_factor + (vt_y * Bz - uz * By)) * path_factor + gx * rel_p * uz
  dpy_ds = (Ey * electric_factor + (uz * Bx - vt_x * Bz)) * path_factor + gy * rel_p * uz

  # Only E changes |p|. dβ/ds = (m c/p₀)² / (E/(p₀ c))³ * d(|p|/p₀)/ds.
  dpz_ds = (Ex * vt_x + Ey * vt_y + Ez * uz) * electric_factor * path_factor
  # dβ/ds * z/β = (m̃² / (rel_p² + m̃²)) * dpz/ds * z/rel_p.
  dz_accel = (tilde_m^2 / kinematics.energy_squared) * dpz_ds * z * inv_rel_p
  dz_ds = beta * reference.inv_beta0 - path_factor + dz_accel

  # Return zero derivatives if momenta are unphysical (branchless)
  zero_deriv = zero(dx_ds)
  return SVector(
    vifelse(good_momenta, dx_ds, zero_deriv),
    vifelse(good_momenta, dpx_ds, zero_deriv),
    vifelse(good_momenta, dy_ds, zero_deriv),
    vifelse(good_momenta, dpy_ds, zero_deriv),
    vifelse(good_momenta, dz_ds, zero_deriv),
    vifelse(good_momenta, dpz_ds, zero_deriv)
  )
end

@inline function _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, gx, gy)
  return _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                      tilde_m, beta_0, gx, gy, _rk_kinematics(pz, tilde_m))
end

@inline function _kick_vector(x, px, y, py, z, pz, s, field::EMField,
                tilde_m, beta_0, gx, gy, kinematics=_rk_kinematics(pz, tilde_m),
                reference=_rk_reference_factors(beta_0, tilde_m))
  Ex, Ey, Ez = field.E
  Bx, By, Bz = field.B
  return _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                     tilde_m, beta_0, gx, gy, kinematics, reference)
end

"""
    kick_vector(x, px, y, py, z, pz, s, field::EMField,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

Evaluate the simplified RK equations with physical electric (V/m) and magnetic
(T) fields. `charge` is in units of e and `p0c` is in eV. The original argument
list is preserved; `mc2` is redundant with `tilde_m` and `p0c` and is unused.
The component overload accepts `Ex, Ey, Ez, Bx, By, Bz` in place of `field`.
"""
@inline function kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                             charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  electric_scale = charge / p0c
  magnetic_scale = electric_scale * c_light(typeof(p0c))
  field = EMField(Ex * magnetic_scale, Ey * magnetic_scale, Ez * magnetic_scale,
                  Bx * magnetic_scale, By * magnetic_scale, Bz * magnetic_scale)
  return _kick_vector(x, px, y, py, z, pz, s, field, tilde_m, beta_0, gx, gy)
end

@inline function kick_vector(x, px, y, py, z, pz, s, field::EMField,
                             charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  return kick_vector(x, px, y, py, z, pz, s, field.E..., field.B...,
                     charge, tilde_m, beta_0, gx, gy, p0c, mc2)
end

"""
  _rk4_step!(coords, i, s, h, field_source, tilde_m, beta_0, gx, gy, magnetic_scale)

Perform a single RK4 step for particle i, updating coordinates in-place.
Only updates state if particle is alive.

# Arguments
- `coords`: Coordinates structure
- `i`: Particle index
- `s`: Current arc length
- `h`: Step size
- `field_source`: Concrete callable field source
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gx`, `gy`: Horizontal and vertical reference curvature components
- `magnetic_scale`: inverse reference rigidity, q/p₀
"""
@inline function _rk4_step!(coords, i, s, h, field_source, tilde_m, beta_0, gx, gy, magnetic_scale,
                            reference=_rk_reference_factors(beta_0, tilde_m))
  # Check if particle is alive
  alive = (coords.state[i] == STATE_ALIVE)
  
  # Extract current particle
  v = coords.v
  x = v[i, XI]
  px = v[i, PXI]
  y = v[i, YI]
  py = v[i, PYI]
  z = v[i, ZI]
  pz = v[i, PZI]

  good = _valid_momentum(px, py, pz)

  # k1 = f(u, s)
  factors = _rk_kinematics(pz, tilde_m)
  t = _rk_time_from_factors(z, s, tilde_m, beta_0, factors, reference)
  field = normalized_field_at(field_source, x, y, s, t, magnetic_scale)
  k1 = _kick_vector(x, px, y, py, z, pz, s, field,
                tilde_m, beta_0, gx, gy, factors, reference)

  # k2 = f(u + h/2 * k1, s + h/2)
  h2 = h / 2
  x2 = x + h2 * k1[1]
  px2 = px + h2 * k1[2]
  y2 = y + h2 * k1[3]
  py2 = py + h2 * k1[4]
  z2 = z + h2 * k1[5]
  pz2 = pz + h2 * k1[6]
  good &= _valid_momentum(px2, py2, pz2)
  factors2 = _rk_kinematics(pz2, tilde_m)
  t2 = _rk_time_from_factors(z2, s + h2, tilde_m, beta_0, factors2, reference)
  field = normalized_field_at(field_source, x2, y2, s + h2, t2, magnetic_scale)
  k2 = _kick_vector(x2, px2, y2, py2, z2, pz2, s + h2, field,
                tilde_m, beta_0, gx, gy, factors2, reference)

  # k3 = f(u + h/2 * k2, s + h/2)
  x3 = x + h2 * k2[1]
  px3 = px + h2 * k2[2]
  y3 = y + h2 * k2[3]
  py3 = py + h2 * k2[4]
  z3 = z + h2 * k2[5]
  pz3 = pz + h2 * k2[6]
  good &= _valid_momentum(px3, py3, pz3)
  factors3 = _rk_kinematics(pz3, tilde_m)
  t3 = _rk_time_from_factors(z3, s + h2, tilde_m, beta_0, factors3, reference)
  field = normalized_field_at(field_source, x3, y3, s + h2, t3, magnetic_scale)
  k3 = _kick_vector(x3, px3, y3, py3, z3, pz3, s + h2, field,
                tilde_m, beta_0, gx, gy, factors3, reference)

  # k4 = f(u + h * k3, s + h)
  x4 = x + h * k3[1]
  px4 = px + h * k3[2]
  y4 = y + h * k3[3]
  py4 = py + h * k3[4]
  z4 = z + h * k3[5]
  pz4 = pz + h * k3[6]
  good &= _valid_momentum(px4, py4, pz4)
  factors4 = _rk_kinematics(pz4, tilde_m)
  t4 = _rk_time_from_factors(z4, s + h, tilde_m, beta_0, factors4, reference)
  field = normalized_field_at(field_source, x4, y4, s + h, t4, magnetic_scale)
  k4 = _kick_vector(x4, px4, y4, py4, z4, pz4, s + h, field,
                tilde_m, beta_0, gx, gy, factors4, reference)

  h6 = h / 6
  xn = x + h6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
  pxn = px + h6 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
  yn = y + h6 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])
  pyn = py + h6 * (k1[4] + 2*k2[4] + 2*k3[4] + k4[4])
  zn = z + h6 * (k1[5] + 2*k2[5] + 2*k3[5] + k4[5])
  pzn = pz + h6 * (k1[6] + 2*k2[6] + 2*k3[6] + k4[6])
  good &= _valid_momentum(pxn, pyn, pzn)
  coords.state[i] = vifelse(alive & !good, STATE_LOST_PZ, coords.state[i])
  accept = alive & good
  v[i, XI] = vifelse(accept, xn, v[i, XI])
  v[i, PXI] = vifelse(accept, pxn, v[i, PXI])
  v[i, YI] = vifelse(accept, yn, v[i, YI])
  v[i, PYI] = vifelse(accept, pyn, v[i, PYI])
  v[i, ZI] = vifelse(accept, zn, v[i, ZI])
  v[i, PZI] = vifelse(accept, pzn, v[i, PZI])
end

"""
    rk4_step!(coords, i, s, h, field_source, charge, tilde_m, beta_0, gx, gy, p0c, mc2)

Advance one RK4 step using a field source with its declared unit convention. The full tracking
kernel reuses the field conversion factors across all steps for each particle.
"""
@inline function rk4_step!(coords, i, s, h, field_source, charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  electric_scale = charge / p0c
  return _rk4_step!(coords, i, s, h, field_source, tilde_m, beta_0, gx, gy,
                   electric_scale * c_light(typeof(p0c)))
end

"""
  rk4_kernel!(i, coords, beta_0, tilde_m, charge, p0c, mc2,
              L, ds_step, n_steps, gx, gy, field_source)

Kernelized RK4 tracking through a concrete electromagnetic field source.
Compatible with @makekernel and the package's kernel architecture.
"""
@makekernel function rk4_kernel!(i, coords::Coords, beta_0, tilde_m, charge, p0c, mc2,
                                L, ds_step, n_steps,
                                gx, gy, field_source)
  s = zero(L)
  electric_scale = charge / p0c
  magnetic_scale = electric_scale * c_light(typeof(p0c))
  reference = _rk_reference_factors(beta_0, tilde_m)

  for step in 1:n_steps
    _rk4_step!(coords, i, s, ds_step, field_source, tilde_m, beta_0, gx, gy, magnetic_scale, reference)
    s += ds_step

    # The common path performs the final callback after exit processing.
    if step != n_steps
      execute_callbacks(i, coords, s, s / (beta_0 * c_light(typeof(ds_step))))
    end
  end
end
