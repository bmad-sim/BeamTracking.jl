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
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, @makekernel, Coords

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

  Ps2 = (1+coords.v[i,PZI])^2 - coords.v[i,PXI]^2 - coords.v[i,PYI]^2        
  coords.state[i] = ifelse(Ps2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i])  

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
  alive = ifelse(coords.state[i] == State.Alive, 1, 0)

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

  v[i,PXI] = alive*(v[i,PXI] * cx - k1 * s * v[i,XI] * sx) - (alive-1) * v[i,PXI]
  v[i,PYI] = alive*(v[i,PYI] * cy + k1 * s * v[i,YI] * sy) - (alive-1) * v[i,PYI]
  v[i,ZI]  = (alive*((v[i,ZI] - (s / 4) * (  xp^2 * (1 + sx * cx)
                                    + yp^2 * (1 + sy * cy)
                                    + k1 / (1 + v[i,PZI])
                                        * ( v[i,XI]^2 * (1 - sx * cx)
                                          - v[i,YI]^2 * (1 - sy * cy) )
                                  )
                      + sign(k1) * ( v[i,XI] * xp * (sqrtks * sx)^2
                              - v[i,YI] * yp * (sqrtks * sy)^2 ) / 2
              )) - (alive-1) * v[i,ZI])
  v[i,XI]  = alive*(v[i,XI] * cx + xp * s * sx) - (alive-1) * v[i,XI]
  v[i,YI]  = alive*(v[i,YI] * cy + yp * s * sy) - (alive-1) * v[i,YI]
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
  PtSqr  = v[i,PXI]^2 + v[i,PYI]^2  # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps2    = P^2 - PtSqr        
  coords.state[i] = ifelse(Ps2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i])
  alive = ifelse(coords.state[i]==State.Alive, 1, 0) 
  Ps = sqrt(Ps2 + (alive-1)*(Ps2-1)) # longitudinal momentum,   Ps   = √[(1 + δ)^2 - P⟂^2]

  v[i,XI] = (alive*(v[i,XI] + s * v[i,PXI] * PtSqr / (P * Ps * (P + Ps)))
            - (alive - 1) * v[i,XI])
  v[i,YI] = (alive*(v[i,YI] + s * v[i,PYI] * PtSqr / (P * Ps * (P + Ps)))
            - (alive - 1) * v[i,YI])
  v[i,ZI] = (alive*(v[i,ZI] - s * ( P * (PtSqr - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(P^2 + tilde_m^2) * Ps
                                    * (beta_0 * sqrt(P^2 + tilde_m^2) + Ps)
                                  )
                            - PtSqr / (2 * P^2)
                          )) - (alive - 1) * v[i,ZI])
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

  Ps2 = (1+coords.v[i,PZI])^2 - coords.v[i,PXI]^2 - coords.v[i,PYI]^2 
  coords.state[i] = ifelse(Ps2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i]) 

  ExactTracking.exact_bend!(      i, coords, e1, e2, g*L/2, g, k0, w, w_inv, tilde_m, beta_0, L / 2)

  if isnothing(coords.q)
    ExactTracking.patch_rotation!(i, coords, w, 0)
    ExactTracking.multipole_kick!(i, coords, mm, knl * 2, ksl * 2, -1)
    ExactTracking.patch_rotation!(i, coords, w_inv, 0)
  else
    ExactTracking.patch_rotation!(i, coords, w, 0)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)
    rotate_spin!(                 i, coords, a, g, tilde_m, mm, kn, ks, L)
    ExactTracking.multipole_kick!(i, coords, mm, knl, ksl, 1)
    ExactTracking.patch_rotation!(i, coords, w_inv, 0)
  end

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
@inline function omega(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  """
  This function computes the integrated spin-precession vector using the multipole 
  coefficients kn and ks indexed by mm, i.e., kn[i] is the normal 
  coefficient of order mm[i].
  """
  @FastGTPSA begin @inbounds begin
    v = coords.v

    # kinetic momenta, not canonical momenta
    px = v[i,PXI] + (v[i,YI] * kn[1] / 2) * (mm[1] == 0) 
    py = v[i,PYI] - (v[i,XI] * kn[1] / 2) * (mm[1] == 0)
  

    rel_p = 1 + v[i,PZI]
    beta_gamma = rel_p / tilde_m
    gamma = sqrt(1 + beta_gamma^2)
    pl2 = rel_p^2 - px^2 - py^2
    coords.state[i] = ifelse(pl2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i])
    alive = ifelse(coords.state[i]==State.Alive, 1, 0) 
    pl = sqrt(pl2 + (alive-1)*(pl2-1))

    beta_hat = SA[px, py, pl] / rel_p

    bx, by = ExactTracking.normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
    bz = kn[1] * (mm[1] == 0)
    b_field = SA[bx, by, bz]

    dot = b_field[1]*beta_hat[1] + b_field[2]*beta_hat[2] + b_field[3]*beta_hat[3]
    b_para = dot * beta_hat
    b_perp = b_field - b_para

    coeff1 = -(1 + g*v[i,XI])/pl
    coeff2 = 1 + a*gamma
    coeff3 = 1 + a

    ox = coeff1 * (coeff2 * b_perp[1] + coeff3 * b_para[1])
    oy = coeff1 * (coeff2 * b_perp[2] + coeff3 * b_para[2]) + g
    oz = coeff1 * (coeff2 * b_perp[3] + coeff3 * b_para[3])

    ox = alive * L * ox
    oy = alive * L * oy
    oz = alive * L * oz

    omega = SA[ox oy oz]
  end end

  #omega = (1 + a*gamma) .* b_perp .+ (1 + a) .* b_para
  #omega = -(1 + g*v[i,XI])/pl .* omega
  #omega = omega .+ SA[0, g, 0]
  #omega = L .* omega
  #omega = alive .* omega
  # The above lines currently give wrong numbers with FastGTPSA

  return omega
end


@makekernel fastgtpsa=true function rotate_spin!(i, coords::Coords, a, g, tilde_m, mm, kn, ks, L)
  """
  This function rotates particle i's quaternion according to the multipoles present.
  """
  q2 = coords.q
  q1 = expq(-omega(i, coords, a, g, tilde_m, mm, kn, ks, L)/2)
  @inbounds begin
    q3 = quat_mul(q1, SA[q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ]])
    q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3[Q0], q3[QX], q3[QY], q3[QZ]
  end
end


@makekernel fastgtpsa=true function integrate_with_spin_thin!(i, coords::Coords, ker, params, a, g, tilde_m, mm, knl, ksl)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  ker(i, coords, params...)
  rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
end

end