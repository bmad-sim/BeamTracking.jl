#=

Tracking using symplectic integration with splits.

=#

macro def_integrator_struct(name)
  quote
    struct $(esc(name))
      order::Int
      num_steps::Int 
      ds_step::Float64
  
      function $(esc(name))(; order::Int=2, num_steps::Int=-1, ds_step::Float64=-1.0)
        _order = order
        _num_steps = num_steps
        _ds_step = ds_step
        if _order ∉ (2, 4, 6, 8)
          error("Symplectic integration only supports orders 2, 4, 6, and 8")
        elseif _num_steps == -1 && _ds_step == -1.0
          _num_steps = 1
        elseif _num_steps > 0 && _ds_step > 0
          error("Only one of num_steps or ds_step should be specified")
        elseif _num_steps < 1 && _ds_step <= 0
          error("Invalid step size")
        elseif _num_steps > 0
          _ds_step = -1.0
        elseif _ds_step > 0
          _num_steps = -1
        end
        return new(_order, _num_steps, _ds_step)
      end
    end
  end
end

@def_integrator_struct(SplitIntegration)
@def_integrator_struct(MatrixKick)
@def_integrator_struct(BendKick)
@def_integrator_struct(SolenoidKick)
@def_integrator_struct(DriftKick)

module IntegrationTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions, ..SIMD, ..SIMDMathFunctions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST, @makekernel, Coords

#
# ===============  I N T E G R A T O R S  ===============
#

@makekernel fastgtpsa=true function order_two_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, L)
  for _ in 1:num_steps
    ker(i, coords, params..., ds_step)
  end
end


@makekernel fastgtpsa=true function order_four_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, L)
  w0 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  w1 =  1.3512071919596577718181151794851757586002349853515625*ds_step
  for _ in 1:num_steps
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1)
  end
end


@makekernel fastgtpsa=true function order_six_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, L)
  w0 =  1.315186320683911169737712043570355*ds_step
  w1 = -1.17767998417887100694641568096432*ds_step
  w2 =  0.235573213359358133684793182978535*ds_step
  w3 =  0.784513610477557263819497633866351*ds_step
  for _ in 1:num_steps
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w3)
  end
end


@makekernel fastgtpsa=true function order_eight_integrator!(i, coords::Coords, ker, params, ds_step, num_steps, L)
  w0 =  1.7084530707869978*ds_step
  w1 =  0.102799849391985*ds_step
  w2 = -1.96061023297549*ds_step
  w3 =  1.93813913762276*ds_step
  w4 = -0.158240635368243*ds_step
  w5 = -1.44485223686048*ds_step
  w6 =  0.253693336566229*ds_step
  w7 =  0.914844246229740*ds_step
  for _ in 1:num_steps
    ker(i, coords, params..., w7)
    ker(i, coords, params..., w6)
    ker(i, coords, params..., w5)
    ker(i, coords, params..., w4)
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w4)
    ker(i, coords, params..., w5)
    ker(i, coords, params..., w6)
    ker(i, coords, params..., w7)
  end
end


#
# ===============  Q U A D R U P O L E  ===============
#
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
@makekernel fastgtpsa=true function mkm_quadrupole!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  rel_p = 1 + coords.v[i,PZI]
  px = coords.v[i,PXI]
  py = coords.v[i,PYI]
  P_s2 = rel_p*rel_p - px*px - py*py
  good_momenta = (P_s2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])

  if !isnothing(coords.q)
    rotate_spin!(               i, coords, a, 0, tilde_m, mm, kn, ks, L / 2)
  end

  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 2)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.patch_rotation!(i, coords, w, 0)
  quadrupole_matrix!(           i, coords, k1, L)
  ExactTracking.patch_rotation!(i, coords, w_inv, 0)
  quadrupole_kick!(             i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 2)

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
         and (signed) reference magnetic R_ref.
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
end # function quadrupole_kick!()


#
# ===============  B E N D  ===============
#

"""
bkb_multipole!()

This integrator uses Bend-Kick-Bend to track a beam through
a curved, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the quadrupole component.

Arguments
—————————
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'k0'       -- dipole strength
- 'mm'       -- order of multipoles
- 'kn'       -- normal multipole strengths 
- 'ks'       -- skew multipole strengths 
- 'L'        -- length
"""
@makekernel fastgtpsa=true function bkb_multipole!(i, coords::Coords, tilde_m, beta_0, a, e1, e2, g, w, w_inv, k0, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  ExactTracking.exact_bend!(      i, coords, e1, e2, g*L/2, g, k0, w, w_inv, tilde_m, beta_0, L / 2)
  ExactTracking.patch_rotation!(  i, coords, w, 0)

  if isnothing(coords.q)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)
    rotate_spin!(                 i, coords, a, g, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)
  end

  ExactTracking.patch_rotation!(  i, coords, w_inv, 0)
  ExactTracking.exact_bend!(      i, coords, e1, e2, g*L/2, g, k0, w, w_inv, tilde_m, beta_0, L / 2)
end 


#
# ===============  S O L E N O I D  ===============
#

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
@makekernel fastgtpsa=true function sks_multipole!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  ExactTracking.exact_solenoid!(i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)

  if isnothing(coords.q)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
    rotate_spin!(                 i, coords, a, 0, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
  end

  ExactTracking.exact_solenoid!(i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
end 


#
# ===============  D R I F T  ===============
#
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
@makekernel fastgtpsa=true function dkd_multipole!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks, L)
  knl = kn * L / 2
  ksl = ks * L / 2

  ExactTracking.exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)

  if isnothing(coords.q)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
  else
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
    rotate_spin!(                 i, coords, a, 0, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, -1)
  end

  ExactTracking.exact_drift!(     i, coords, beta_0, gamsqr_0, tilde_m, L / 2)
end


#
# ===============  S P I N  ===============
#

# WARNING!!! IF YOU INLINE THIS FUNCTION, SPIN TRACKING THROUGH A SOLENOID WILL
# BREAK WITH SIMD!!!
# BEWARE!!!
"""
This function computes the integrated spin-precession vector using the multipole 
coefficients kn and ks indexed by mm, i.e., knl[i] is the normal 
coefficient of order mm[i].
"""
function omega(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  @FastGTPSA begin
    v = coords.v

    # kinetic momenta, not canonical momenta
    if (length(kn) == 0) || (mm[1] != 0)
      px = v[i,PXI]
      py = v[i,PYI] 
    else
      px = v[i,PXI] + (v[i,YI] * kn[1] / 2)
      py = v[i,PYI] - (v[i,XI] * kn[1] / 2)
    end

    rel_p = 1 + v[i,PZI]
    beta_gamma = rel_p / tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    pl2 = rel_p*rel_p - px*px - py*py
    pl2_0 = zero(pl2)
    good_momenta = (pl2 > pl2_0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    pl2_1 = one(pl2)
    pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

    bx, by = ExactTracking.normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
    bz = mm[1] == 0 ? kn[1] : zero(kn[1])

    coeff = -(1 + g*v[i,XI])/pl
    coeff1 = coeff * (1 + a*gamma)
    coeff2 = coeff * (1 + a)

    dot = bx*px/rel_p + by*py/rel_p + bz*pl/rel_p

    b_para_x = dot * px / rel_p
    b_para_y = dot * py / rel_p
    b_para_z = dot * pl / rel_p

    b_perp_x = (bx - b_para_x) * coeff1
    b_perp_y = (by - b_para_y) * coeff1
    b_perp_z = (bz - b_para_z) * coeff1

    b_para_x = b_para_x * coeff2
    b_para_y = b_para_y * coeff2
    b_para_z = b_para_z * coeff2

    ox = (b_perp_x + b_para_x) * L        
    oy = (b_perp_y + b_para_y + g) * L
    oz = (b_perp_z + b_para_z) * L

    omega = (ox, oy, oz)
  end

  return omega
end


"""
This function rotates particle i's quaternion according to the multipoles present.
"""
@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega(i, coords, a, g, tilde_m, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end


@makekernel fastgtpsa=true function integrate_with_spin_thin!(i, coords::Coords, ker, params, a, g, tilde_m, mm, knl, ksl)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  ker(i, coords, params...)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
end

end