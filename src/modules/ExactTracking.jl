#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix
# (equal to number of temporaries needed for a single particle)
struct Exact end

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
using ..BeamTracking: C_LIGHT
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@makekernel fastgtpsa=true function update_P0!(i, b, Brho_initial, Brho_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * Brho_initial / Brho_final
  end
  return v
end
=#

#
# ===============  E X A C T   D R I F T  ===============
#
"""
exact_drift!()

In the computation of z_final, we use the fact that
    1/√a - 1/√b == (b - a)/(√a √b (√a + √b))
to avoid the potential for severe cancellation when
a and b both have the form 1 + ε for different small
values of ε.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
L: element length
"""
@makekernel fastgtpsa=true function exact_drift!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
  v = b.v
  P_s = sqrt((1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2)) 
  v[i,XI]   = v[i,XI] + v[i,PXI] * L / P_s
  v[i,YI]   = v[i,YI] + v[i,PYI] * L / P_s
  # high-precision computation of z_final
  # vf.z = vi.z - (1 + δ) * L * (1 / Ps - 1 / (β0 * sqrt((1 + δ)^2 + tilde_m^2)))
  v[i,ZI]   = v[i,ZI] - ( (1 + v[i,PZI]) * L
                * ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                / ( beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) * P_s
                    * (beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) + P_s)
                  )
              )
end # function exact_drift!()


#
# ===============  Q U A D R U P O L E  ===============
#
"""
mkm_quadrupole!()

This integrator uses Matrix-Kick-Matrix to implement a quadrupole
integrator accurate though second-order in the integration step-size.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
k2_num:   g / Bρ0 = g / (p0 / q)
          where g and Bρ0 respectively denote the quadrupole gradient
          and (signed) reference magnetic rigidity.
L: element length
"""
@makekernel fastgtpsa=true function mkm_quadrupole!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, k2_num, L)
  quadrupole_matrix!(i, b::BunchView, k2_num, L / 2)
  quadrupole_kick!(  i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
  quadrupole_matrix!(i, b::BunchView, k2_num, L / 2)
end # function mkm_quadrupole!()


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k2_num:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic rigidity.
s: element length
"""
@makekernel fastgtpsa=true function quadrupole_matrix!(i, b::BunchView, k2_num, s)
  v = b.v

  sgn = sign(k2_num)
  focus = k2_num >= 0  # horizontally focusing for positive particles

  xp = v[i,PXI] / (1 + v[i,PZI])  # x'
  yp = v[i,PYI] / (1 + v[i,PZI])  # y'
  sqrtks = sqrt(abs(k2_num / (1 + v[i,PZI]))) * s  # |κ|s
  cx = focus ? cos(sqrtks) : cosh(sqrtks)
  cy = focus ? cosh(sqrtks) : cos(sqrtks)
  sx = focus ? sincu(sqrtks) : sinhcu(sqrtks)
  sy = focus ? sinhcu(sqrtks) : sincu(sqrtks)

  v[i,PXI] = v[i,PXI] * cx - k2_num * s * v[i,XI] * sx
  v[i,PYI] = v[i,PYI] * cy + k2_num * s * v[i,YI] * sy
  v[i,ZI]  = (v[i,ZI] - (s / 4) * (  xp^2 * (1 + sx * cx)
                                    + yp^2 * (1 + sy * cy)
                                    + k2_num / (1 + v[i,PZI])
                                        * ( v[i,XI]^2 * (1 - sx * cx)
                                          - v[i,YI]^2 * (1 - sy * cy) )
                                  )
                      + sgn * ( v[i,XI] * xp * (sqrtks * sx)^2
                              - v[i,YI] * yp * (sqrtks * sy)^2 ) / 2
              )
  v[i,XI]  = v[i,XI] * cx + xp * s * sx
  v[i,YI]  = v[i,YI] * cy + yp * s * sy
end # function quadrupole_matrix!()


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
  rP0 = 1 + v[i,PZI]              # reduced total momentum,  P/P0 = 1 + δ
  sqrPt = v[i,PXI]^2 + v[i,PYI]^2   # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  Ps = sqrt(rP0^2 - sqrPt)          # longitudinal momentum,   Ps = √[(1 + δ)^2 - P⟂^2]
  v[i,XI] = v[i,XI] + s * v[i,PXI] / rP0 * sqrPt / (Ps * (rP0 + Ps))
  v[i,YI] = v[i,YI] + s * v[i,PYI] / rP0 * sqrPt / (Ps * (rP0 + Ps))
  v[i,ZI] = v[i,ZI] - s * ( rP0
                              * (sqrPt - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(rP0^2 + tilde_m^2) * Ps
                                    * (beta_0 * sqrt(rP0^2 + tilde_m^2) + Ps)
                                  )
                            - sqrPt / (2 * rP0^2)
                          )
end # function quadrupole_kick!()


#
# ===============  M U L T I P O L E  ===============
#
"""
dkd_multipole()

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
  exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
  multipole_kick!(i, b, mm, kn * L, ks * L)
  exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L / 2)
end # function dkd_multipole!()


"""
    multipole_kick!(i, b, ms, knl, ksl)

Track a beam of particles through a thin-lens multipole
having integrated normal and skew strengths listed in the
coefficient vectors knl and ksl respectively. The vector ms
lists the orders of the corresponding entries in knl and ksl.

The algorithm used in this function takes advantage of the
complex representation of the vector potential Az,
  - ``-Re{ sum_m (b_m + i a_m) (x + i y)^m / m }``,
and uses a Horner-like scheme (see Shachinger and Talman
[SSC-52]) to compute the transverse kicks induced by a pure
multipole magnet. This method supposedly has good numerical
properties, though I've not seen a proof of that claim.

DTA: Ordering matters!
DTA: Add thin dipole kick.

### Arguments
 - ms:  vector of m values for non-zero multipole coefficients
 - knl: vector of normal integrated multipole strengths
 - ksl: vector of skew integrated multipole strengths


     NB: Here the j-th component of knl (ksl) denotes the
       normal (skew) component of the multipole strength of
       order mm[j] (after scaling by the reference Bρ).
       For example, if mm[j] = 3, then knl[j] denotes the
       normal integrated sextupole strength scaled by Bρo.
"""
#@inline function multipole_kick!(i, b, mm, knl, ksl)
#  @inbounds begin #@FastGTPSA! begin
#  jm = length(mm)
#  m = mm[jm]
#  if m == 1 return v end
#  work[i,1] = ar = knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
#  work[i,2] = ai = knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
#  jm -= 1
#  while m > 2
#    m -= 1
#    work[i,3] = work[i,1] * v[i,XI] - work[i,2] * v[i,YI]
#    work[i,2] = work[i,1] * v[i,YI] + work[i,2] * v[i,XI]
#    work[i,1] = work[i,3]
#    if m == mm[jm]
#      work[i,1] += knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
#      work[i,2] += knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
#      jm -= 1
#    end
#  end
#  v[i,PXI] -= work[i,1]
#  v[i,PYI] += work[i,2]
#  end #end
#  return v
#end # function multipole_kick!()
#
@makekernel fastgtpsa=true function multipole_kick!(i, b::BunchView, ms, knl, ksl)
  v = b.v
  jm = length(ms)
  m  = ms[jm]
  knl_tot = knl[jm]
  ksl_tot = ksl[jm]
  jm -= 1
  while 2 <= m
    m -= 1
    tmp_knl_tot = knl_tot * v[i,XI] - ksl_tot * v[i,YI]
    ksl_tot     = knl_tot * v[i,YI] + ksl_tot * v[i,XI]
    knl_tot = tmp_knl_tot
    if 0 < jm && m == ms[jm]
      knl_tot += knl[jm]
      ksl_tot += ksl[jm]
      jm -= 1
    end
  end
  v[i,PXI] -= knl_tot
  v[i,PYI] += ksl_tot
end # function multipole_kick!()


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
"""
This function implements exact symplectic tracking through a
sector bend, derived using the Hamiltonian (25.9) given in the
BMad manual. As a consequence of using that Hamiltonian, the
reference value of βγ must be that of a particle with the
design energy.  Should we wish to change that, we shall need
to carry both reference and design values.

Arguments
—————————
beta_0: β_0 = (βγ)_0 / √(γ_0^2)
brho_0: Bρ_0,  reference magnetic rigidity
hc: coordinate frame curvature
b0: magnet field strength
e1: entrance face angle (+ve angle <=> toward rbend)
e2: exit face angle (+ve angle <=> toward rbend)
Lr: element arc length
"""
@makekernel fastgtpsa=true function exact_sbend!(i, b::BunchView, beta_0, brho_0, hc, b0, e1, e2, Lr)
  v = b.v

  rho = brho0 / b0
  ang = hc * Lr
  c1 = cos(ang)
  s1 = sin(ang)

  P_s     = sqrt((1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
  P_alpha = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2)                 # P_α
  s1phx = (1 + hc * v[i,XI]) / (hc * rho)                     # scaled (1 + h x)
  Pxpph = P_s - s1phx                                 # Px'/h
  ang_eff = ang + asin(v[i,PXI] / P_alpha) - asin((v[i,PXI] * c1 + Pxpph * s1) / P_alpha)  # α + φ1 - φ2
  # high-precision computation of x-final
  v[i,XI] = (v[i,XI] * c1 - Lr * sin(ang / 2) * sincu(ang / 2)
             + rho * (v[i,PXI] + ((v[i,PXI]^2 + (P_s + Pxpph) * s1phx) * s1 - 2v[i,PXI] * Pxpph * c1)
                             / (sqrt(P_alpha^2 - (v[i,PXI] * c1 + Pxpph * s1)^2) + P_s * c1)) * s1)
  v[i,PXI] = v[i,PXI] * c1 + Pxpph * s1
  v[i,YI] = v[i,YI] + rho * v.py * ang_eff
  # high-precision computation of z-final
  v[i,ZI] = (v[i,ZI] - rho * (1 + v[i,PZI]) * ang_eff
               + (1 + v[i,PZI]) * Lr / (beta_0 * sqrt(1 / beta_0^2 + (2 + v[i,PZI]) * v[i,PZI])))
end # function exact_sbend!()


@makekernel fastgtpsa=true function exact_solenoid!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, L)
  v = b.v
  # Recurring variables 
  rel_p = 1 + v[i,PZI]    
  pr = sqrt(rel_p^2 - (v[i,PXI] + v[i,YI] * ks / 2)^2 - (v[i,PYI] - v[i,XI] * ks / 2)^2)
  s = sin(ks * L / pr)   
  cp = 1 + cos(ks * L / pr)                
  cm = 2 - cp
  # Temporaries
  x_0 = v[i,XI] 
  px_0 = v[i,PXI] 
  y_0 = v[i,YI]  
  # Update
  v[i,ZI]  -= rel_p * L *
                ((v[i,PXI] + v[i,YI] * ks / 2)^2 + (v[i,PYI] - v[i,XI] * ks / 2)^2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) /
                ( beta_0 * sqrt(rel_p^2 + tilde_m^2) * pr * (beta_0 * sqrt(rel_p^2 + tilde_m^2) + pr) )
  v[i,XI] = cp * x_0 / 2 + s * (px_0 / ks + y_0 / 2) + cm * v[i,PYI] / ks
  v[i,PXI] = s * (v[i,PYI] / 2 - ks * x_0 / 4) + cp * px_0 / 2 - ks * cm * y_0 / 4
  v[i,YI] = s * (v[i,PYI] / ks - x_0 / 2) + cp * y_0 / 2 - cm * px_0 / ks
  v[i,PYI]  = ks * cm * x_0 / 4 - s * (px_0 / 2 + ks * y_0 / 4) + cp * v[i,PYI] / 2
end

@makekernel fastgtpsa=true function patch!(i, b::BunchView, tilde_m, dt, dx, dy, dz, winv::Union{StaticMatrix{3,3},Nothing}, L)
  # Temporary momentum [1+δp, ps_0]
  v = b.v
  rel_p = 1 + v[i,PZI]                                 
  ps_0 = sqrt(rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2)  
  # Only apply rotations if needed
  if isnothing(winv)
    # No rotation case
    v[i,XI] -= dx
    v[i,YI] -= dy

    # Apply t_offset
    v[i,ZI] += rel_p/sqrt(rel_p^2+tilde_m^2)*C_LIGHT*dt

    # Drift to face
    v[i,XI]   += v[i,PXI] * dz / ps_0
    v[i,YI]   += v[i,PYI] * dz / ps_0
    v[i,ZI]   -=  dz * rel_p / ps_0 - L*rel_p*sqrt((1+tilde_m^2)/(rel_p^2+tilde_m^2))
  else
    # Translate position vector [x, y]
    x_0 = v[i,XI] - dx                                # x_0
    y_0 = v[i,YI] - dy                                # y_0

    # Temporary momentum vector [px, py]
    px_0 = v[i,PXI]                                    # px_0
    py_0 = v[i,PYI]                                    # py_0

    # Transform position vector [x - dx, y - dy, -dz]
    v[i,XI]   = winv[1,1]*x_0 + winv[1,2]*y_0 - winv[1,3]*dz
    v[i,YI]   = winv[2,1]*x_0 + winv[2,2]*y_0 - winv[2,3]*dz
    s_f = winv[3,1]*x_0 + winv[3,2]*y_0 - winv[3,3]*dz  # s_f

    # Transform momentum vector [px, py, ps]
    v[i,PXI]  = winv[1,1]*px_0 + winv[1,2]*py_0 + winv[1,3]*ps_0
    v[i,PYI]  = winv[2,1]*px_0 + winv[2,2]*py_0 + winv[2,3]*ps_0
    ps_f = winv[3,1]*px_0 + winv[3,2]*py_0 + winv[3,3]*ps_0 # ps_f

    # Apply t_offset
    v[i,ZI] += rel_p/sqrt(rel_p^2+tilde_m^2)*C_LIGHT*dt

    # Drift to face
    v[i,XI] -= s_f * v[i,PXI] / ps_f
    v[i,YI] -= s_f * v[i,PYI] / ps_f
    v[i,ZI] += s_f * rel_p / ps_f + L*rel_p*sqrt((1+tilde_m^2)/(rel_p^2+tilde_m^2))
  end
end


# Utility functions ============================================================

# Rotation matrix
"""
  w_matrix(x_rot, y_rot, z_rot)

Constructs a rotation matrix based on the given Bryan-Tait angles.

Bmad/SciBmad follows the MAD convention of applying z, x, y rotations in that order.
Furthermore, in ReferenceFrameRotations, the rotation angles are defined as negative
of the SciBmad rotation angles `x_rot`, `y_rot`, and `z_rot`.

The inverse matrix reverses the order of operations and their signs.


Arguments:
- `x_rot::Number`: Rotation angle around the x-axis.
- `y_rot::Number`: Rotation angle around the y-axis.
- `z_rot::Number`: Rotation angle around the z-axis.

Returns:
- `DCM{Float64}`: ReferenceFrameRotations.DCM (direct cosine matrix), rotation matrix.
"""
function w_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(-z_rot, -x_rot, -y_rot, :ZXY)
end

# Inverse rotation matrix
function w_inv_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(y_rot, x_rot, z_rot, :YXZ)
end

function drift_params(species::Species, Brho)
  beta_gamma_0 = BeamTracking.calc_beta_gammma(species, Brho)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0^2
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end


"""
    exact_bend!(i, v, work, theta, gtot, g, L, mc2, p0c, beta_ref) 

Tracks a particle through a sector bend via exact tracking. (no edge angles)

#Arguments
- 'theta'    -- 'g' * 'L'
- 'gtot'     -- 'g' + 'dg'
- 'p0c'      -- reference momentum in eV
- 'beta_ref' -- 'p0c' / sqrt('mc2'^2 + 'p0c'^2)
"""
@inline function exact_bend!(i, v, work, theta, gtot, g, L, mc2, p0c, beta_ref)
  @inbounds begin #@FastGTPSA! begin
      work[i,1] = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2) #pt
      work[i,2] = theta + asin(v[i,PXI] / work[i,1]) #phi1
      work[i,3] = gtot / work[i,1] #gp
      work[i,4] = 1+g*v[i,XI] # 1 + g * x
      work[i,5] = cos(work[i,2]) #cos(theta + phi1)
      work[i,6] = sincu(theta) #sincu(theta)
      work[i,7] = 2*work[i,4]*sin(work[i,2])*L*work[i,6]- work[i,3]*work[i,4]^2*L^2*(work[i,6])^2 #alpha
      if abs(work[i,2]) < π/2
          work[i,8] = work[i,7]/(sqrt(work[i,5]^2 + work[i,3]*work[i,7]) + work[i,5]) #xi
      else
          work[i,8] = (sqrt(work[i,5]^2 + work[i,3]*work[i,7]) - work[i,5]) / work[i,3] #xi
      end
      work[i,9] = -L*work[i,6]-v[i,XI]*sin(theta) #Lcv
      work[i,10] = 2 * (work[i,2] - atan(work[i,8], -work[i,9])) #theta_p
      work[i,11] = sqrt(work[i,9]^2 + work[i,8]^2) / sincu(work[i,10]/2) #Lp
      work[i,12] = (v[i,PZI] * p0c + p0c) #p

      v[i,XI] = v[i,XI]*cos(theta) - g/2*L^2*(sincu(theta/2))^2 + work[i,8]
      v[i,PXI] = work[i,1]*sin(work[i,2] - work[i,10])
      v[i,YI] = v[i,YI] + v[i,PYI]*work[i,11]/work[i,1] 
      v[i,ZI] = v[i,ZI] - (1 + v[i,PZI])*work[i,11]/work[i,1] + 
                      L*(work[i,12]/sqrt(mc2^2+work[i,12]^2))/beta_ref
  end #end
  return v
end

end
