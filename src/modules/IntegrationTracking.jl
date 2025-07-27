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
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView

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
@makekernel fastgtpsa=true function mkm_quadrupole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, w, w_inv, k1, mm, kn, ks, L)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)
  quadrupole_kick!(             i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.patch_rotation!(i, b, w, 0)
  quadrupole_matrix!(           i, b, k1, L)
  ExactTracking.patch_rotation!(i, b, w_inv, 0)
  quadrupole_kick!(             i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)
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
@makekernel fastgtpsa=true function sks_multipole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, Ksol, mm, kn, sn, L)
  ExactTracking.exact_solenoid!(i, b, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, sn * L, 1)
  ExactTracking.exact_solenoid!(i, b, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
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
starting with the dipole component. (For example, b[3] denotes
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
@makekernel fastgtpsa=true function dkd_multipole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, mm, kn, ks, L)
  ExactTracking.exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L, ks * L, 1)
  ExactTracking.exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
end


#
# ===============  IMPLICIT  ===============
#

@makekernel fastgtpsa=true function tilted_kik_multipole!(i, b::BunchView, β0, tilde_m, w, w_inv, gx, gy, A, mm, kn, ks, L)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
  ExactTracking.patch_rotation!(i, b, w, 0)
  implicit_step!(               i, b, gx, gy, A, L)
  ExactTracking.patch_rotation!(i, b, w_inv, 0)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 3)
end 

@makekernel fastgtpsa=true function kik_multipole!(i, b::BunchView, β0, tilde_m, gx, gy, A, mm, kn, ks, L)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
  implicit_step!(               i, b, gx, gy, A, L)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
  ExactTracking.multipole_kick!(i, b, mm, kn * L / 2, ks * L / 2, 1)
end 

@makekernel fastgtpsa=true function did_implicit!(i, b::BunchView, β0, tilde_m, gx, gy, A, L)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
  implicit_step!(               i, b, gx, gy, A, L)
  implicit_drift!(              i, b, tilde_m, β0, L / 2)
end 


@makekernel fastgtpsa=true function implicit_step!(i, b::BunchView, gx, gy, A, ds)
  @inline implicit_step_1!(i, b, gx, gy, A, ds / 2)
  @inline implicit_step_2!(i, b, gx, gy, A, ds / 2)
end

@makekernel fastgtpsa=true function implicit_drift!(i, b::BunchView, tilde_m, β0, ds)
  v = b.v
  rel_p = 1 + v[i,PZI]

  v[i,ZI] += ds * rel_p / (β0 * sqrt(rel_p^2 + tilde_m^2))
end

"""
implicit_step_1!()

Arguments
—————————
gx: horizontal curvature
gy: vertical curvature
A: (x,y) -> (ax, ay) transverse vector potential function (longitudinally constant)
ds: step size

"""
function implicit_step_1!(i, b::BunchView, gx, gy, A, ds)
  v = b.v

  T = eltype(v)
  xdual = (
    BeamTracking.Dual(v[i,XI], (one(T), zero(T))),
    BeamTracking.Dual(v[i,YI], (zero(T), one(T))),
  )
  ydual = A(xdual)
  ax, ay = BeamTracking.value.(ydual)
  jac = (
    BeamTracking.partials.(ydual, 1),
    BeamTracking.partials.(ydual, 2),
  )

  function dF(zf; z0) 
    px = zf[PXI] - ax
    py = zf[PYI] - ay
    rel_p = 1 + zf[PZI]
    ps = sqrt(rel_p^2 - px^2 - py^2)
    dps = ds * (1 + gx * z0[XI] + gy * z0[YI]) / ps
    dF = SVector{6}(
      zf[XI]  - z0[XI]  - dps * px,
      zf[YI]  - z0[YI]  - dps * py,
      zf[ZI]  - z0[ZI]  + dps * rel_p,
      zf[PXI] - z0[PXI] + dps * (jac[1][1]*px + jac[2][1]*py) + gx * ps * ds,
      zf[PYI] - z0[PYI] + dps * (jac[1][2]*px + jac[2][2]*py) + gy * ps * ds,
      zf[PZI] - z0[PZI]
    )
    return dF
  end
  BeamTracking.newton!(dF, @view v[i,:])
end

"""
implicit_step!()

Arguments
—————————
gx: horizontal curvature
gy: vertical curvature
A: (x,y) -> (ax, ay) transverse vector potential function (longitudinally constant)
ds: step size

"""
function implicit_step_2!(i, b::BunchView, gx, gy, A, ds)
  v = b.v

  function dF(zf; z0) 
    T = eltype(zf)
    xdual = (
      BeamTracking.Dual(zf[XI], (one(T), zero(T))),
      BeamTracking.Dual(zf[YI], (zero(T), one(T)))
    )
    ydual = A(xdual)
    ax, ay = BeamTracking.value.(ydual)
    jac = (
      BeamTracking.partials.(ydual, 1),
      BeamTracking.partials.(ydual, 2)
    )

    px = z0[PXI] - ax
    py = z0[PYI] - ay
    rel_p = 1 + z0[PZI]
    ps = sqrt(rel_p^2 - px^2 - py^2)
    dps = ds * (1 + gx * zf[XI] + gy * zf[YI]) / ps
    dF = SVector{6}(
      zf[XI]  - z0[XI]  - dps * px,
      zf[YI]  - z0[YI]  - dps * py,
      zf[ZI]  - z0[ZI]  + dps * rel_p,
      zf[PXI] - z0[PXI] + dps * (jac[1][1]*px + jac[2][1]*py) + gx * ps * ds,
      zf[PYI] - z0[PYI] + dps * (jac[1][2]*px + jac[2][2]*py) + gy * ps * ds,
      zf[PZI] - z0[PZI]
    )
    return dF
  end
  BeamTracking.newton!(dF, @view v[i,:])
end

function implicit_full_1!(i, b::BunchView, A, tilde_m, β0, ds)
  v = b.v

  function H(z)
    a = A(z[[XI,YI]])
    return -sqrt((1+z[PZI])^2 - (z[PXI]-a[1])^2 - (z[PYI]-a[2])^2) - a[3] + sqrt((1+z[PZI])^2 + tilde_m^2) / β0
  end
  
  function dF(zf; z0) 
    T = eltype(zf)
    xdual = (
      BeamTracking.Dual(z0[XI], ( one(T), zero(T), zero(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(zf[PXI],(zero(T),  one(T), zero(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(z0[YI], (zero(T), zero(T),  one(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(zf[PYI],(zero(T), zero(T), zero(T),  one(T), zero(T), zero(T))),
      BeamTracking.Dual(z0[ZI], (zero(T), zero(T), zero(T), zero(T),  one(T), zero(T))),
      BeamTracking.Dual(zf[PZI],(zero(T), zero(T), zero(T), zero(T), zero(T),  one(T))),
    )
    ydual = H(xdual)
    H_diff = (
      BeamTracking.partials.(ydual, 1),
      BeamTracking.partials.(ydual, 2),
      BeamTracking.partials.(ydual, 3),
      BeamTracking.partials.(ydual, 4),
      BeamTracking.partials.(ydual, 5),
      BeamTracking.partials.(ydual, 6),
    )

    dF = SVector{6}(
      zf[XI]  - z0[XI]  - ds * H_diff[2],
      zf[YI]  - z0[YI]  - ds * H_diff[4],
      zf[ZI]  - z0[ZI]  - ds * H_diff[6],
      zf[PXI] - z0[PXI] + ds * H_diff[1],
      zf[PYI] - z0[PYI] + ds * H_diff[3],
      zf[PZI] - z0[PZI] + ds * H_diff[5]
    )
    return dF
  end
  BeamTracking.newton!(dF, @view v[i,:])
end

function implicit_full_2!(i, b::BunchView, A, tilde_m, β0, ds)
  v = b.v

  function H(z)
    a = A(z[[XI,YI]])
    return -sqrt((1+z[PZI])^2 - (z[PXI]-a[1])^2 - (z[PYI]-a[2])^2) - a[3] + sqrt((1+z[PZI])^2 + tilde_m^2) / β0
  end

  function dF(zf; z0) 
    T = eltype(zf)
    xdual = (
      BeamTracking.Dual(zf[XI], ( one(T), zero(T), zero(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(z0[PXI],(zero(T),  one(T), zero(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(zf[YI], (zero(T), zero(T),  one(T), zero(T), zero(T), zero(T))),
      BeamTracking.Dual(z0[PYI],(zero(T), zero(T), zero(T),  one(T), zero(T), zero(T))),
      BeamTracking.Dual(zf[ZI], (zero(T), zero(T), zero(T), zero(T),  one(T), zero(T))),
      BeamTracking.Dual(z0[PZI],(zero(T), zero(T), zero(T), zero(T), zero(T),  one(T))),
    )
    ydual = H(xdual)
    H_diff = (
      BeamTracking.partials.(ydual, 1),
      BeamTracking.partials.(ydual, 2),
      BeamTracking.partials.(ydual, 3),
      BeamTracking.partials.(ydual, 4),
      BeamTracking.partials.(ydual, 5),
      BeamTracking.partials.(ydual, 6),
    )

    dF = SVector{6}(
      zf[XI]  - z0[XI]  - ds * H_diff[2],
      zf[YI]  - z0[YI]  - ds * H_diff[4],
      zf[ZI]  - z0[ZI]  - ds * H_diff[6],
      zf[PXI] - z0[PXI] + ds * H_diff[1],
      zf[PYI] - z0[PYI] + ds * H_diff[3],
      zf[PZI] - z0[PZI] + ds * H_diff[5]
    )
    return dF
  end
  BeamTracking.newton!(dF, @view v[i,:])
end

@makekernel fastgtpsa=true function implicit_full_step!(i, b::BunchView, A, tilde_m, β0, ds)
  implicit_full_1!(i, b, A, tilde_m, β0, ds/2)
  implicit_full_2!(i, b, A, tilde_m, β0, ds/2)
end


function H(z; A, tilde_m, β0)
  @FastGTPSA begin
    a = A(z[[XI,YI]])
    return -sqrt((1+z[PZI])^2 - (z[PXI]-a[1])^2 - (z[PYI]-a[2])^2) - a[3] + sqrt((1+z[PZI])^2 + tilde_m^2) / β0
  end
end

function F(z::Vector{TPS64{D}}; A, tilde_m=1.0, β0=1.0, ds=1, p::Bool=true) where D
  @FastGTPSA begin
  h = H(z; A=A, tilde_m=tilde_m, β0=β0)
  return z + (2*p-1) * ds * [
          GTPSA.deriv(h,2);
          GTPSA.deriv(h,1);
          GTPSA.deriv(h,4);
          GTPSA.deriv(h,3);
          GTPSA.deriv(h,6);
          GTPSA.deriv(h,5)]
  end
end

@makekernel fastgtpsa=true function implicit_TPSA_step_1!(i, b::BunchView, A, tilde_m, β0, ds)
  f = BunchView([copy(b.state)[i]], GTPSA.scalar.(b.v[i,:]'), nothing)
  implicit_full_2!(i, f, A, tilde_m, β0, ds)
  
  D = eltype(b.v).parameters[2]
  coord = GTPSA.scalar.([f.v[1, XI], 
                         b.v[i, PXI], 
                         f.v[1, YI], 
                         b.v[i, PYI], 
                         f.v[1, ZI], 
                         b.v[i, PZI]]) + @vars(D)
  f1 = F(coord; A=A, tilde_m=tilde_m, β0=β0, ds=ds, p=false)

  f1_pinv = zero(f1)  # Partial inversion
  GTPSA.pminv!(6, f1, 6, f1_pinv, [1, 0, 1, 0, 1, 0])

  b.v[i,:] = f.v[1,:] + f1_pinv ∘ cutord.(b.v[i,:], 0)  # Compose GTPS and add scalar part
end

@makekernel fastgtpsa=true function implicit_TPSA_step_2!(i, b::BunchView, A, tilde_m, β0, ds)
  f = BunchView([copy(b.state)[i]], GTPSA.scalar.(b.v[i,:]'), nothing)
  implicit_full_1!(i, f, A, tilde_m, β0, ds)

  D = eltype(b.v).parameters[2] # GTPSA descriptor
  coord = GTPSA.scalar.([ b.v[i, XI], 
                          f.v[1, PXI], 
                          b.v[i, YI], 
                          f.v[1, PYI], 
                          b.v[i, ZI], 
                          f.v[1, PZI]]) + @vars(D)
  f2 = F(coord; A=A, tilde_m=tilde_m, β0=β0, ds=ds, p=true)

  f2_pinv = zero(f2)  # Partial inversion
  GTPSA.pminv!(6, f2, 6, f2_pinv, [0, 1, 0, 1, 0, 1])

  b.v[i,:] = f.v[1,:] + f2_pinv ∘ cutord.(b.v[i,:], 0)  # Compose GTPS and add scalar part
end

@makekernel fastgtpsa=true function implicit_TPSA_step!(i, b::BunchView, A, tilde_m, β0, ds)
  implicit_TPSA_step_1!(i, b, A, tilde_m, β0, ds/2)
  implicit_TPSA_step_2!(i, b, A, tilde_m, β0, ds/2)
end

end