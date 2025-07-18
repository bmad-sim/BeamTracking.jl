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
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, @makekernel, BunchView

#
# ===============  I N T E G R A T O R S  ===============
#

@makekernel fastgtpsa=true function order_two_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  for _ in 1:num_steps
    ker(i, b, params..., ds_step)
  end
end


@makekernel fastgtpsa=true function order_four_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  w0 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  w1 = 1.3512071919596577718181151794851757586002349853515625*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
  end
end


@makekernel fastgtpsa=true function order_six_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  w0 = 1.315186320683911169737712043570355*ds_step
  w1 = -1.17767998417887100694641568096432*ds_step
  w2 = 0.235573213359358133684793182978535*ds_step
  w3 = 0.784513610477557263819497633866351*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
  end
end


@makekernel fastgtpsa=true function order_eight_integrator!(i, b::BunchView, ker, params, ds_step, num_steps, L)
  w0 = 1.7084530707869978*ds_step
  w1 = 0.102799849391985*ds_step
  w2 = -1.96061023297549*ds_step
  w3 = 1.93813913762276*ds_step
  w4 = -0.158240635368243*ds_step
  w5 = -1.44485223686048*ds_step
  w6 = 0.253693336566229*ds_step
  w7 = 0.914844246229740*ds_step
  for _ in 1:num_steps
    ker(i, b, params..., w7)
    ker(i, b, params..., w6)
    ker(i, b, params..., w5)
    ker(i, b, params..., w4)
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
    ker(i, b, params..., w4)
    ker(i, b, params..., w5)
    ker(i, b, params..., w6)
    ker(i, b, params..., w7)
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
@makekernel fastgtpsa=true function mkm_quadrupole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, beta_gamma_0, a, w, w_inv, k1, mm, kn, ks, L)
  if !isnothing(b.q)
    rotate_spin!(               i, b, a, 0, beta_gamma_0, mm, kn, ks, L / 2)
  end

  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)
  quadrupole_kick!(             i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.patch!(         i, b, tilde_m, 0, 0, 0, 0, w, 0)
  quadrupole_matrix!(           i, b, k1, L)
  ExactTracking.patch!(         i, b, tilde_m, 0, 0, 0, 0, w_inv, 0)
  quadrupole_kick!(             i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)

  if !isnothing(b.q)
    rotate_spin!(               i, b, a, 0, beta_gamma_0, mm, kn, ks, L / 2)
  end
end 


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k1:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic rigidity.
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_matrix!(i, b::BunchView, k1, s)
  v = b.v

  focus = k1 >= 0  # horizontally focusing if positive
  defocus = k1 < 0 

  rel_p = 1 + v[i,PZI]
  xp = v[i,PXI] / rel_p  # x'
  yp = v[i,PYI] / rel_p  # y'
  sqrtks = sqrt(abs(k1 / rel_p)) * s  # |κ|s

  cosine = cos(sqrtks) # precompute trig for branchless but still fast
  coshine = cosh(sqrtks)
  sinecu = sincu(sqrtks)
  shinecu = sinhcu(sqrtks)
  cx = focus * cosine  + defocus * coshine # branchless
  cy = focus * coshine + defocus * cosine
  sx = focus * sinecu  + defocus * shinecu
  sy = focus * shinecu + defocus * sinecu

  v[i,PXI] = v[i,PXI] * cx - k1 * s * v[i,XI] * sx
  v[i,PYI] = v[i,PYI] * cy + k1 * s * v[i,YI] * sy
  v[i,ZI]  = (v[i,ZI] - (s / 4) * (  xp^2 * (1 + sx * cx)
                                    + yp^2 * (1 + sy * cy)
                                    + k1 / (1 + v[i,PZI])
                                        * ( v[i,XI]^2 * (1 - sx * cx)
                                          - v[i,YI]^2 * (1 - sy * cy) )
                                  )
                      + sign(k1) * ( v[i,XI] * xp * (sqrtks * sx)^2
                              - v[i,YI] * yp * (sqrtks * sy)^2 ) / 2
              )
  v[i,XI]  = v[i,XI] * cx + xp * s * sx
  v[i,YI]  = v[i,YI] * cy + yp * s * sy
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
@makekernel fastgtpsa=true function quadrupole_kick!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, s)
  v = b.v

  P     = 1 + v[i,PZI]             # [scaled] total momentum, P/P0 = 1 + δ
  PtSqr = v[i,PXI]^2 + v[i,PYI]^2  # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps    = sqrt(P^2 - PtSqr)        # longitudinal momentum,   Ps   = √[(1 + δ)^2 - P⟂^2]

  v[i,XI] = v[i,XI] + s * v[i,PXI] * PtSqr / (P * Ps * (P + Ps))
  v[i,YI] = v[i,YI] + s * v[i,PYI] * PtSqr / (P * Ps * (P + Ps))
  v[i,ZI] = v[i,ZI] - s * ( P * (PtSqr - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(P^2 + tilde_m^2) * Ps
                                    * (beta_0 * sqrt(P^2 + tilde_m^2) + Ps)
                                  )
                            - PtSqr / (2 * P^2)
                          )
end # function quadrupole_kick!()


#=
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
beta_0: β_0 = (βγ)_0 / √(γ_0^2)
brho_0: Bρ_0,  reference magnetic rigidity
hc: coordinate frame curvature
b0: dipole field strength
e1: entrance face angle (+ve angle <=> toward rbend)
e2: exit face angle (+ve angle <=> toward rbend)
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
sn: vector of skew multipole strengths scaled by Bρ0
L:  element arc length
"""
@makekernel fastgtpsa=true function bkb_multipole!(i, b::BunchView, beta_0, brho_0, hc, b0, e1, e2, mm, kn, sn, L)
  ExactTracking.exact_sbend!(   i, b, beta_0, brho_0, hc, b0, e1, e2, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, sn * L)
  ExactTracking.exact_sbend!(   i, b, beta_0, brho_0, hc, b0, e1, e2, L / 2)
end 
=#


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
@makekernel fastgtpsa=true function sks_multipole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, beta_gamma_0, a, Ksol, mm, kn, ks, L)
  ExactTracking.exact_solenoid!(  i, b, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)

  if isnothing(b.q)
    ExactTracking.multipole_kick!(i, b, mm, kn * L, ks * L, 1)
  else
    ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
    rotate_spin!(                 i, b, a, 0, beta_gamma_0, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
  end

  ExactTracking.exact_solenoid!(  i, b, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
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
starting with the dipole component.

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
@makekernel fastgtpsa=true function dkd_multipole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, beta_gamma_0, a, mm, kn, ks, L)
  ExactTracking.exact_drift!(     i, b, beta_0, gamsqr_0, tilde_m, L / 2)

  if isnothing(b.q)
    ExactTracking.multipole_kick!(i, b, mm, kn * L, ks * L, 1)
  else
    ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
    rotate_spin!(                 i, b, a, 0, beta_gamma_0, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
  end

  ExactTracking.exact_drift!(     i, b, beta_0, gamsqr_0, tilde_m, L / 2)
end


#
# ===============  S P I N  ===============
#
@inline function omega(i, b::BunchView, a, g, beta_gamma_0, mm, kn, ks)
  """
  This function computes the spin-precession vector using the multipole 
  coefficients kn and ks indexed by mm, i.e., knl[i] is the normal 
  coefficient of order mm[i].
  """
  v = b.v

  px = v[i,PXI] + (v[i,YI] * kn[1] / 2) * (mm[1] == 0) # kinetic momenta
  py = v[i,PYI] - (v[i,XI] * kn[1] / 2) * (mm[1] == 0)

  rel_p = 1 + v[i,PZI]
  gamma_0 = sqrt(1 + beta_gamma_0^2)
  beta_gamma = rel_p * beta_gamma_0
  gamma = sqrt(1 + beta_gamma^2)
  pl = sqrt(rel_p^2 - px^2 - py^2)
  beta_hat = SA[px, py, pl] / rel_p

  bx, by = ExactTracking.normalized_field!(mm, kn, ks, v[i,XI], v[i,YI], 1)
  bz = kn[1] * (mm[1] == 0)
  b_field = SA[bx, by, bz]

  dot = b_field[1]*beta_hat[1] + b_field[2]*beta_hat[2] + b_field[3]*beta_hat[3]
  b_para = dot * beta_hat
  b_perp = b_field - b_para
  
  omega = (1 + a*gamma)*b_perp + (1 + a)*b_para
  omega = -(1 + g*v[i,XI])/pl * omega
  omega = omega + SA[0, g, 0]

  return omega
end


@inline function expq(v)
  """
  This function computes exp(i v⋅σ) as a quaternion, where σ is the 
  vector of Pauli matrices.
  """
  n = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
  c = cos(n)
  s = sincu(n)
  v2 = s * v
  return SA[-c, v2[1], v2[2], v2[3]]
end


@makekernel fastgtpsa=true function rotate_spin!(i, b::BunchView, a, g, beta_gamma_0, mm, kn, ks, L)
  """
  This function rotates b.q according to the multipoles present.
  """
  q2 = b.q
  q1 = expq(-L / 2 * omega(i, b, a, g, beta_gamma_0, mm, kn, ks))
  a1, b1, c1, d1 = q1[1], q1[2], q1[3], q1[4]
  a2, b2, c2, d2 = q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ]
  q2[i,Q0] = a1*a2 - b1*b2 - c1*c2 - d1*d2
  q2[i,QX] = a1*b2 + b1*a2 + c1*d2 - d1*c2
  q2[i,QY] = a1*c2 - b1*d2 + c1*a2 + d1*b2
  q2[i,QZ] = a1*d2 + b1*c2 - c1*b2 + d1*a2
end


@makekernel fastgtpsa=true function integrate_with_spin_thin!(i, b::BunchView, ker, params, a, g, beta_gamma_0, mm, knl, ksl)
  rotate_spin!(i, b, a, g, beta_gamma_0, mm, knl, ksl)
  ker(i, b, params...)
end


end