#=

  “Exact” tracking methods

=#

struct Exact end

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions, ..SIMDMathFunctions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST, @makekernel, Coords, vifelse, BeamTracking.track_rotation!
using ..BeamTracking: C_LIGHT
const TRACKING_METHOD = Exact


@makekernel fastgtpsa=true function update_P0!(i, coords::Coords, R_ref_initial, R_ref_final, ramp_without_rf)
  v = coords.v 

  v[i,PXI] = v[i,PXI] * R_ref_initial / R_ref_final
  v[i,PYI] = v[i,PYI] * R_ref_initial / R_ref_final
  if !ramp_without_rf
    v[i,PZI] = (R_ref_initial * (1 + v[i,PZI]) - R_ref_final) / R_ref_final
  end
end


#
# ===============  E X A C T   D R I F T  ===============
#
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
end # function exact_drift!()


#
# ===============  M U L T I P O L E  ===============
#
"""
    multipole_kick!(i, coords, ms, knl, ksl)

Track a beam of particles through a thin-lens multipole
having integrated normal and skew strengths listed in the
coefficient vectors knl and ksl respectively. The vector ms
lists the orders of the corresponding entries in knl and ksl.

The algorithm used in this function takes advantage of the
complex representation of the vector potential Az,
  - ``-Re{ sum_m (b_m + i a_m) (x + i y)^m / m! }``,
and uses a Horner-like scheme (see Shachinger and Talman
[SSC-52]) to compute the transverse kicks induced by a pure
multipole magnet. This method supposedly has good numerical
properties, though I've not seen a proof of that claim.

## Arguments
 - ms:  vector of m values for non-zero multipole coefficients
 - knl: vector of normal integrated multipole strengths
 - ksl: vector of skew integrated multipole strengths


     NB: Here the j-th component of knl (ksl) denotes the
       normal (skew) component of the multipole strength of
       order ms[j] (after scaling by the reference Bρ).
       For example, if ms[j] = 3, then knl[j] denotes the
       normal integrated sextupole strength scaled by Bρo.
       Moreover, and this is essential, the multipole
       coefficients must appear in ascending order.
"""
@makekernel fastgtpsa=true function multipole_kick!(i, coords::Coords, ms, knl, ksl, excluding)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  bx, by = normalized_field(ms, knl, ksl, v[i,XI], v[i,YI], excluding)
  bx_0 = zero(bx)
  by_0 = zero(by)
  v[i,PXI] -= vifelse(alive, by, by_0)                   
  v[i,PYI] += vifelse(alive, bx, bx_0)
end # function multipole_kick!()


function normalized_field(ms, knl, ksl, x, y, excluding)
  """
  Returns (bx, by), the transverse components of the magnetic field divided
  by the reference rigidty.
  """
  @FastGTPSA begin
    jm = length(ms)
    m  = ms[jm]
    add = (m != excluding && m > 0)
    knl_0 = zero(knl[jm]*x)
    ksl_0 = zero(ksl[jm]*y)
    by_0 = knl[jm]*one(x)
    bx_0 = ksl[jm]*one(y)
    by = vifelse(add, by_0, knl_0)
    bx = vifelse(add, bx_0, ksl_0)
    jm -= 1
    while 2 <= m
      m -= 1
      t  = (by * x - bx * y) / m
      bx = (by * y + bx * x) / m
      by = t
      add = (0 < jm && m == ms[jm]) && (m != excluding) # branchless
      idx = max(1, jm) # branchless trickery
      new_by = by + knl[idx]
      new_bx = bx + ksl[idx]
      by = vifelse(add, new_by, by)
      bx = vifelse(add, new_bx, bx)
      jm -= add
    end
  end
  return bx, by
end


#function binom(m::Integer, x, y)
#  """
#  This function computes the real and imaginary parts of
#  (x + i y)^m. One can use these components to compute,
#  among other things, the multipole kick induced by two-
#  dimensional multipole magnets.
#  """
#  if m == 0
#    return [ 1, 0.0 ]
#  end
#  ar = x
#  ai = y
#  mm = m
#  while mm > 1
#    mm -= 1
#    t  = x * ar - y * ai
#    ai = y * ar + x * ai
#    ar = t
#  end
#  return [ ar, ai ]
#end # function binom()


#
# ===============  E X A C T   S E C T O R   B E N D  ===============
#
#=
"""
    exact_sbend!(i, coords, β0, Bρ0, hc, b0, e1, e2, Larc)
This function implements exact symplectic tracking through a
sector bend, derived using the Hamiltonian (25.9) given in the
BMad manual. As a consequence of using that Hamiltonian, the
reference value of βγ must be that of a particle with the
design energy.  Should we wish to change that, we shall need
to carry both reference and design values.

## Arguments
- beta_0: β_0 = (βγ)_0 / √(γ_0^2)
- R_ref:  Bρ_0,  reference magnetic R_ref
- hc:     coordinate frame curvature
- b0:     magnet field strength
- e1:     entrance face angle (+ve angle <=> toward rbend)
- e2:     exit face angle (+ve angle <=> toward rbend)
- Larc:   element arc length, in meters
"""
@makekernel fastgtpsa=true function exact_sbend!(i, coords::Coords, beta_0, R_ref, hc, b0, e1, e2, Lr)
  v = coords.v

  rho = R_ref0 / b0
  ang = hc * Lr
  c1 = cos(ang)
  s1 = sin(ang)

  P_s     = sqrt((1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
  P_alpha = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2)                 # P_α
  s1phx = (1 + hc * v[i,XI]) / (hc * rho)                     # scaled (1 + h x)
  Pxpph = P_s - s1phx                                 # Px'/h
  ang_eff = ang + asin(v[i,PXI] / P_alpha) - asin((v[i,PXI] * c1 + Pxpph * s1) / P_alpha)  # α + φ1 - φ2

  # high-precision computation of x-final:
  v[i,XI] = (v[i,XI] * c1 - Lr * sin(ang / 2) * sincu(ang / 2)
             + rho * (v[i,PXI] + ((v[i,PXI]^2 + (P_s + Pxpph) * s1phx) * s1 - 2v[i,PXI] * Pxpph * c1)
                             / (sqrt(P_alpha^2 - (v[i,PXI] * c1 + Pxpph * s1)^2) + P_s * c1)) * s1)
  v[i,PXI] = v[i,PXI] * c1 + Pxpph * s1
  v[i,YI] = v[i,YI] + rho * v.py * ang_eff

  # high-precision computation of z-final
  v[i,ZI] = (v[i,ZI] - rho * (1 + v[i,PZI]) * ang_eff
               + (1 + v[i,PZI]) * Lr / (beta_0 * sqrt(1 / beta_0^2 + (2 + v[i,PZI]) * v[i,PZI])))
end # function exact_sbend!()
=#

"""
    exact_bend!(i, coords::Coords, e1, e2, theta, g, Kn0, w, w_inv, tilde_m, beta_0, L)

Tracks a particle through a sector bend via exact tracking. If edge angles are 
provided, a linear hard-edge fringe map is applied at both ends.

#Arguments
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'Kn0'      -- normalized dipole field
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'L'        -- length
"""
@makekernel fastgtpsa=true function exact_bend!(i, coords::Coords, theta, g, Kn0, tilde_m, beta_0, L)
  v = coords.v
  rel_p = 1 + v[i,PZI]

  pt2 = rel_p*rel_p - v[i,PYI]*v[i,PYI]
  good_momenta = (pt2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pt2_1 = one(pt2)
  pt = sqrt(vifelse(good_momenta, pt2, pt2_1))

  arg = v[i,PXI] / pt
  abs_arg = abs(arg)
  arg_1 = one(arg)
  arg_0 = zero(arg)
  good_arg = (abs_arg < arg_1)
  coords.state[i] = vifelse(!good_arg & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)

  phi1 = theta + asin(vifelse(good_arg, arg, arg_0))
  gp = Kn0 / pt
  h = 1 + g*v[i,XI] 
  cplus = cos(phi1) 
  splus = sin(phi1)
  sinc_theta = sincu(theta)
  sinc_theta_2 = sincu(theta/2)
  cosc_theta = sinc_theta_2*sinc_theta_2/2
  sgn = sign(L)
  alpha_helper = h*L*sinc_theta
  alpha = 2*h*splus*L*sinc_theta - gp*alpha_helper*alpha_helper

  cond = cplus*cplus + gp*alpha
  good_cond = (cond > 0)
  coords.state[i] = vifelse(!good_cond & alive, STATE_LOST, coords.state[i]) # particle does not intersect the exit face
  alive = (coords.state[i] == STATE_ALIVE)
  cond_1 = one(cond)
  nasty_sqrt = sqrt(vifelse(good_cond, cond, cond_1))

  gp_0 = zero(gp)
  abs_gp = abs(gp)
  good_gp = (abs_gp > gp_0)
  gp_1 = one(gp)
  gp_safe = vifelse(good_gp, gp, gp_1)
  pos_cplus = (cplus > 0)
  xi1 = alpha/(nasty_sqrt + cplus)
  xi2 = (nasty_sqrt - cplus)/gp_safe
  xi = vifelse(!good_gp | pos_cplus, xi1, xi2)

  Lcv = -sgn*(L*sinc_theta + v[i,XI]*sin(theta)) 
  negative_Lcv = -Lcv
  thetap = 2*(phi1 - sgn*atan2(xi, negative_Lcv)) 
  Lp = sgn*sqrt(Lcv*Lcv + xi*xi)/sincu(thetap/2) 

  new_x = v[i,XI]*cos(theta) - L*L*g*cosc_theta + xi
  new_px = pt*sin(phi1 - thetap)
  new_y = v[i,YI] + v[i,PYI]*Lp/pt
  new_z = v[i,ZI] - rel_p*Lp/pt + L*rel_p/sqrt(tilde_m*tilde_m + rel_p*rel_p)/beta_0
  v[i,XI]  = vifelse(alive, new_x, v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y, v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z, v[i,ZI])
end


@makekernel fastgtpsa=true function linear_bend_fringe!(i, coords::Coords, f)
  v = coords.v

  alive = (coords.state[i] == STATE_ALIVE)
  new_px = v[i,PXI] + v[i,XI]f
  new_py = v[i,PYI] - v[i,YI]f
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
end



@makekernel fastgtpsa=true function exact_bend_with_rotation!(i, coords::Coords, e1, e2, theta, g, Kn0, w, w_inv, tilde_m, beta_0, L)
  BeamTracking.track_rotation!(i, coords, w, 0)
  ExactTracking.linear_bend_fringe!(i, coords, Kn0*tan(e1))
  exact_bend!(i, coords, theta, g, Kn0, tilde_m, beta_0, L)
  ExactTracking.linear_bend_fringe!(i, coords, Kn0*tan(e2))
  BeamTracking.track_rotation!(i, coords, w_inv, 0)
end


# This is separate because the spin can be transported exactly here
@makekernel fastgtpsa=true function exact_curved_drift!(i, coords::Coords, e1, e2, theta, g, w, w_inv, a, tilde_m, beta_0, L) 
  exact_bend_with_rotation!(i, coords, 0, 0, theta, g, 0, w, w_inv, tilde_m, beta_0, L)
  if !isnothing(coords.q)
    BeamTracking.track_rotation!(i, coords, w, 0)
    IntegrationTracking.rotate_spin!(i, coords, a, g, tilde_m, SA[0], SA[0], SA[0], L)
    BeamTracking.track_rotation!(i, coords, w_inv, 0)
  end
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


@makekernel fastgtpsa=true function patch_offset!(i, coords::Coords, tilde_m, dx, dy, dt)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p = 1 + v[i,PZI]
  new_x = v[i,XI] - dx
  new_y = v[i,YI] - dy
  new_z = v[i,ZI] + rel_p/sqrt(rel_p*rel_p+tilde_m*tilde_m)*C_LIGHT*dt
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end

@makekernel fastgtpsa=true function patch!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, winv, L) 
  v = coords.v
  rel_p = 1 + v[i,PZI]
  ps_02 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  ps_02_0 = zero(ps_02)
  good_momenta = (ps_02 > ps_02_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_02_1 = one(ps_02)
  ps_0 = sqrt(vifelse(good_momenta, ps_02, ps_02_1))
  # Only apply rotations if needed
  if isnothing(winv)
    patch_offset!(i, coords, tilde_m, dx, dy, dt)
    exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, L)
    new_z = v[i,ZI] - (dz - L) * rel_p / ps_0
    v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
  else
    patch_offset!(i, coords, tilde_m, dx, dy, dt)
    w31 = 2*(winv[QX]*winv[QZ] - winv[QY]*winv[Q0])
    w32 = 2*(winv[QY]*winv[QZ] + winv[QX]*winv[Q0])
    w33 = 1 - 2*(winv[QX]*winv[QX] + winv[QY]*winv[QY])
    s_f = w31*v[i,XI] + w32*v[i,YI] - w33*dz
    BeamTracking.track_rotation!(i, coords, winv, -dz)
    exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, -s_f)
    new_z = v[i,ZI] + ((s_f + L) * rel_p * 
    sqrt((1 + tilde_m*tilde_m)/(rel_p*rel_p + tilde_m*tilde_m)))
    v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
  end
end


# Utility functions ============================================================

function drift_params(species::Species, R_ref)
  beta_gamma_0 = BeamTracking.R_to_beta_gamma(species, R_ref)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0*beta_gamma_0
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end

end