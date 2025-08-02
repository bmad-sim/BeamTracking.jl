#=

  “Exact” tracking methods

=#

struct Exact end

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
using ..BeamTracking: C_LIGHT
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@makekernel fastgtpsa=true function update_P0!(i, b, rigidity_initial, rigidity_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * rigidity_initial / rigidity_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * rigidity_initial / rigidity_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * rigidity_initial / rigidity_final
  end
  return v
end
=#

#
# ===============  E X A C T   D R I F T  ===============
#
"""
    exact_drift!(i, b, β_0, γsqr_0, tilde_m, L)

Return the result of exact tracking a particle through a drift
of length `L`, assuming `β_0`, `γsqr_0`, and `tilde_m` respectively
denote the reference velocity normalized to the speed of light,
the corresponding value of the squared Lorentz factor, and the
particle rest energy normalized to the reference value of ``pc``.

NB: In the computation of ``z_final``, we use the fact that
  - ``1/√a - 1/√b == (b - a)/(√a √b (√a + √b))``
to avoid the potential for severe cancellation when
``a`` and ``b`` both have the form ``1 + ε`` for different small
values of ``ε``.

## Arguments
- `β_0`:     reference velocity normalized to the speed of light, ``v_0 / c``
- `γsqr_0`:  corresponding value of the squared Lorentz factor
- `tilde_m`: particle rest energy normalized to the reference value of ``pc``
- `L`:       element length, in meters
"""
@makekernel fastgtpsa=true function exact_drift!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
  v = b.v

  P_s2 = (1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2)
  b.state[i] = ifelse(P_s2 <= 0 && b.state[i] == State.Alive, State.Lost, b.state[i])
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 
  P_s = sqrt(P_s2 + (alive-1)*(P_s2-1))

  v[i,XI] = alive*(v[i,XI] + v[i,PXI] * L / P_s) - (alive - 1) * v[i,XI]
  v[i,YI] = alive*(v[i,YI] + v[i,PYI] * L / P_s) - (alive - 1) * v[i,YI]
  # high-precision computation of z_final:
  #   vf.z = vi.z - (1 + δ) * L * (1 / Ps - 1 / (β0 * sqrt((1 + δ)^2 + tilde_m^2)))
  v[i,ZI] = alive*(v[i,ZI] - ( (1 + v[i,PZI]) * L
                * ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                / ( beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) * P_s
                    * (beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) + P_s)
                  )
              )) - (alive - 1) * v[i,ZI]
end # function exact_drift!()


#
# ===============  M U L T I P O L E  ===============
#
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
@makekernel fastgtpsa=true function multipole_kick!(i, b::BunchView, ms, knl, ksl, start)
  v = b.v

  jm = length(ms)
  m  = ms[jm]
  add = (start <= m)
  ar = knl[jm] * add
  ai = ksl[jm] * add
  jm -= 1
  while 2 <= m
    m -= 1
    t  = (ar * v[i,XI] - ai * v[i,YI]) / m
    ai = (ar * v[i,YI] + ai * v[i,XI]) / m
    ar = t
    add = (0 < jm && m == ms[jm]) && (start <= m) # branchless
    idx = max(1, jm) # branchless trickery
    ar += knl[idx] * add
    ai += ksl[idx] * add
    jm -= add
  end
  v[i,PXI] -= ar
  v[i,PYI] += ai
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
    exact_sbend!(i, b, β0, Bρ0, hc, b0, e1, e2, Larc)
This function implements exact symplectic tracking through a
sector bend, derived using the Hamiltonian (25.9) given in the
BMad manual. As a consequence of using that Hamiltonian, the
reference value of βγ must be that of a particle with the
design energy.  Should we wish to change that, we shall need
to carry both reference and design values.

## Arguments
- beta_0: β_0 = (βγ)_0 / √(γ_0^2)
- rigidity_0: Bρ_0,  reference magnetic rigidity
- hc:     coordinate frame curvature
- b0:     magnet field strength
- e1:     entrance face angle (+ve angle <=> toward rbend)
- e2:     exit face angle (+ve angle <=> toward rbend)
- Larc:   element arc length, in meters
"""
@makekernel fastgtpsa=true function exact_sbend!(i, b::BunchView, beta_0, rigidity_0, hc, b0, e1, e2, Lr)
  v = b.v

  rho = rigidity0 / b0
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

"""
    exact_bend!(i, b::BunchView, e1, e2, theta, g, Kn0, w, w_inv, tilde_m, beta_0, L)

Tracks a particle through a sector bend via exact tracking. If edge angles are 
provided, a linear hard-edge fringe map is applied at both ends.

#Arguments
- 'e1'       -- entrance face angle
- 'e2'       -- exit face angle
- 'theta'    -- 'g' * 'L'
- 'g'        -- curvature
- 'Kn0'      -- normalized dipole field
- 'w'        -- rotation matrix into curvature/field plane
- 'w_inv'    -- rotation matrix out of curvature/field plane
- 'tilde_m'  -- mc2/p0c
- 'beta_0'   -- p0c/E0
- 'L'        -- length
"""
@makekernel fastgtpsa=false function exact_bend!(i, b::BunchView, e1, e2, theta, g, Kn0, w::StaticMatrix{3,3}, w_inv::StaticMatrix{3,3}, tilde_m, beta_0, L)
  me1 = Kn0*tan(e1)/(1 + b.v[i,PZI])
  mx1 = SA[1 0; me1  1]
  my1 = SA[1 0;-me1  1]
  me2 = Kn0*tan(e2)/(1 + b.v[i,PZI])
  mx2 = SA[1 0; me2  1]
  my2 = SA[1 0;-me2 1]
  
  patch_rotation!(i, b, w, 0)
  LinearTracking.linear_coast_uncoupled!(i, b, mx1, my1, 0, nothing, nothing)

  v = b.v
  rel_p = 1 + v[i,PZI]
 
  pt2 = rel_p^2 - v[i,PYI]^2
  b.state[i] = ifelse(pt2 <= 0 && b.state[i] == State.Alive, State.Lost, b.state[i])
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 
  
  pt = sqrt(pt2 + (alive-1)*(pt2-1))
  arg = v[i,PXI] / pt
  b.state[i] = ifelse(abs(arg) > 1 && b.state[i] == State.Alive, State.Lost, b.state[i])
  # The above comparison does not work with FastGTPSA (currently)
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 

  phi1 = theta + asin(alive * arg)
  gp = Kn0 / pt
  h = 1 + g*v[i,XI] 
  cplus = cos(phi1) 
  splus = sin(phi1)
  sinc_theta = sincu(theta)
  cosc_theta = (sincu(theta/2))^2 / 2
  sgn = sign(L)
  alpha = 2*h*splus*L*sinc_theta - gp*(h*L*sinc_theta)^2

  cond = cplus^2 + gp*alpha
  b.state[i] = ifelse(cond <= 0 && b.state[i] == State.Alive, State.Lost, b.state[i]) # particle does not intersect the exit face
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 
  nasty_sqrt = alive * sqrt(cond + (alive-1)*(cond-1))

  xi = ifelse(cplus > 0 || gp ≈ 0, alpha/(nasty_sqrt + cplus), (nasty_sqrt - cplus)/(gp + ((abs(gp)>0)-1)*(gp-1)))

  Lcv = -sgn*(L*sinc_theta + v[i,XI]*sin(theta)) 
  thetap = 2 * (phi1 - sgn*atan(xi, -Lcv)) 
  Lp = sgn*sqrt(Lcv^2 + xi^2)/sincu(thetap/2) 

  v[i,XI]  = alive*(v[i,XI]*cos(theta) - L^2*g*cosc_theta + xi) - (alive - 1) * v[i,XI]
  v[i,PXI] = alive*(pt*sin(phi1 - thetap)) - (alive - 1) * v[i,PXI]
  v[i,YI]  = alive*(v[i,YI] + v[i,PYI]*Lp/pt) - (alive - 1) * v[i,YI]
  v[i,ZI]  = alive*(v[i,ZI] - rel_p*Lp/pt + 
                  L*rel_p/sqrt(tilde_m^2+rel_p^2)/beta_0) - (alive - 1) * v[i,ZI]

  LinearTracking.linear_coast_uncoupled!(i, b, mx2, my2, 0, nothing, nothing)
  patch_rotation!(i, b, w_inv, 0)
end

@makekernel fastgtpsa=true function exact_solenoid!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, L)
  v = b.v

  # Recurring variables
  rel_p = 1 + v[i,PZI]
  pr2 = rel_p^2 - (v[i,PXI] + v[i,YI] * ks / 2)^2 - (v[i,PYI] - v[i,XI] * ks / 2)^2
  b.state[i] = ifelse(pr2 <= 0 && b.state[i] == State.Alive, State.Lost, b.state[i])
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 
  pr = sqrt(pr2 + (alive-1)*(pr2-1))

  s = sin(ks * L / pr)
  cp = 1 + cos(ks * L / pr)
  cm = 2 - cp
  # Temporaries
  x_0 = v[i,XI]
  px_0 = v[i,PXI]
  y_0 = v[i,YI]
  
  # Update
  v[i,ZI]  = (alive*(v[i,ZI] - rel_p * L *
                ((v[i,PXI] + v[i,YI] * ks / 2)^2 + (v[i,PYI] - v[i,XI] * ks / 2)^2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) /
                ( beta_0 * sqrt(rel_p^2 + tilde_m^2) * pr * (beta_0 * sqrt(rel_p^2 + tilde_m^2) + pr) )) 
                - (alive - 1) * v[i,ZI])
  v[i,XI] = alive*(cp * x_0 / 2 + s * (px_0 / ks + y_0 / 2) + cm * v[i,PYI] / ks) - (alive - 1) * v[i,XI]
  v[i,PXI] = alive*(s * (v[i,PYI] / 2 - ks * x_0 / 4) + cp * px_0 / 2 - ks * cm * y_0 / 4) - (alive - 1) * v[i,PXI]
  v[i,YI] = alive*(s * (v[i,PYI] / ks - x_0 / 2) + cp * y_0 / 2 - cm * px_0 / ks) - (alive - 1) * v[i,YI]
  v[i,PYI]  = alive*(ks * cm * x_0 / 4 - s * (px_0 / 2 + ks * y_0 / 4) + cp * v[i,PYI] / 2) - (alive - 1) * v[i,PYI]
end

@makekernel fastgtpsa=true function patch_offset!(i, b::BunchView, tilde_m, dx, dy, dt)
  v = b.v
  rel_p = 1 + v[i,PZI]
  v[i,XI] -= dx
  v[i,YI] -= dy
  v[i,ZI] += rel_p/sqrt(rel_p^2+tilde_m^2)*C_LIGHT*dt
end

@makekernel fastgtpsa=true function patch_rotation!(i, b::BunchView, winv::StaticMatrix{3,3}, dz)
  v = b.v
  cond = (1 + v[i,PZI])^2 - v[i,PXI]^2 - v[i,PYI]^2
  b.state[i] = ifelse(cond <= 0 && b.state[i] == State.Alive, State.Lost, b.state[i])
  alive = ifelse(b.state[i]==State.Alive, 1, 0) 
  ps_0 = alive * sqrt(cond + (alive-1)*(cond-1))
  x_0 = v[i,XI]
  y_0 = v[i,YI]
  v[i,XI]   = alive*(winv[1,1]*x_0 + winv[1,2]*y_0 - winv[1,3]*dz) - (alive - 1)*v[i,XI]
  v[i,YI]   = alive*(winv[2,1]*x_0 + winv[2,2]*y_0 - winv[2,3]*dz) - (alive - 1)*v[i,YI]

  px_0 = v[i,PXI]
  py_0 = v[i,PYI]
  v[i,PXI] = alive*(winv[1,1]*px_0 + winv[1,2]*py_0 + winv[1,3]*ps_0) - (alive - 1)*v[i,PXI]
  v[i,PYI] = alive*(winv[2,1]*px_0 + winv[2,2]*py_0 + winv[2,3]*ps_0) - (alive - 1)*v[i,PYI]
end

@makekernel fastgtpsa=true function patch!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, winv::Union{StaticMatrix{3,3},Nothing}, L)
  v = b.v
  rel_p = 1 + v[i,PZI]
  # Only apply rotations if needed
  if isnothing(winv)
    patch_offset!(i, b, tilde_m, dx, dy, dt)
    exact_drift!(i, b, beta_0, gamsqr_0, tilde_m, L)
    v[i,ZI] -= (dz-L) * rel_p / sqrt(rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2)
  else
    patch_offset!(i, b, tilde_m, dx, dy, dt)
    s_f = winv[3,1]*v[i,XI] + winv[3,2]*v[i,YI] - winv[3,3]*dz
    patch_rotation!(i, b, winv, dz)
    exact_drift!(i, b, beta_0, gamsqr_0, tilde_m, -s_f)
    v[i,ZI] += (s_f + L) * rel_p * sqrt((1 + tilde_m^2)/(rel_p^2 + tilde_m^2))
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

function drift_params(species::Species, rigidity)
  beta_gamma_0 = BeamTracking.calc_beta_gamma(species, rigidity)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0^2
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end

end
