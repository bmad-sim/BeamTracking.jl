module MatrixKick
using ..GTPSA: @FastGTPSA!, GTPSA
using StructArrays
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

Base.@kwdef struct Drift{T}
  L::T  # drift length / m
end

Base.@kwdef struct Quadrupole{T}
  L::T    # quadrupole length / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end


function track!(bunch::Bunch, ele::MatrixKick.Drift; work=get_work(bunch, Val{1}()))
#=
This function implements symplectic tracking through a drift,
derived using the Hamiltonian (25.9) given in the BMad manual.
As a consequence of using that Hamiltonian, the reference value
of βγ must be that of a particle with the design energy.
Should we wish to change that, we shall need to carry both
reference and design values.
=#
  L = ele.L
  v = bunch.v

  tilde_m    = 1 / bunch.beta_gamma_ref
  gamsqr_ref = 1 + bunch.beta_gamma_ref^2
  beta_ref   = bunch.beta_gamma_ref / sqrt(gamsqr_ref)

  @FastGTPSA! begin
  @. work[1] = sqrt((1.0 + v.pz)^2 - (v.px^2 + v.py^2))  # P_s
  @. v.x  .= v.x + v.px * L / work[1]
  @. v.y  .= v.y + v.py * L / work[1]
  #@. v.z  .= v.z - ( (1.0 + v.pz) * L
  #                   * (1. / work[1] - 1. / (beta_ref * sqrt((1.0 + v.pz)^2 + tilde_m^2))) )
  # high-precision computation of z-final
  @. v.z  .= v.z - ( (1.0 + v.pz) * L
                     * ((v.px^2 + v.py^2) - v.pz * (2 + v.pz) / gamsqr_ref)
                     / ( beta_ref * sqrt((1.0 + v.pz)^2 + tilde_m^2) * work[1]
                         * (beta_ref * sqrt((1.0 + v.pz)^2 + tilde_m^2) + work[1])
                       )
                   )
  end

  # Spin unchanged

  return bunch
end # function track!(::Bunch, ::Drift)



# This integrator uses the so-called Matrix-Kick-Matrix method to implement
# an integrator accurate though second-order in the integration step-size.
function track!(bunch::Bunch, ele::MatrixKick.Quadrupole; work=get_work(bunch, Val{6}()))
  L = ele.L

  # κ^2 (kappa-squared) := (q g / P0) / (1 + δ)
  # numerator of κ^2
  k2_num = ele.Bn1 / brho(massof(bunch.species), bunch.beta_gamma_ref, chargeof(bunch.species))

  v = bunch.v
  v_work = StructArray{Coord{eltype(work[1])}}((work[1], work[2], work[3], work[4], work[5], work[6]))

  trackQuadMx!(v_work, v, k2_num, L / 2)
  trackQuadK!( v, v_work, bunch.beta_gamma_ref, L)
  trackQuadMx!(v_work, v, k2_num, L / 2)

  v .= v_work
  return bunch
end # function track!(::Bunch, ::Quadrupole)


"""
    trackQuadMx!(vf::StructArray, vi::StructArray, k2_num::Float64, s::Float64)

track "matrix part" of quadrupole
"""
function trackQuadMx!(vf, vi, k2_num, s)
  focus   = k2_num >= 0  # horizontally focusing
  defocus = k2_num <  0  # horizontally defocusing

  p =  @. 1 + vi.pz    # reduced momentum, P/P0 = 1 + δ
  k2 = @. k2_num / p  # κ^2 for each particle
  ks = @. sqrt(abs(k2)) * s  # |κ|s
  xp = @. vi.px / p  # x'
  yp = @. vi.py / p  # y'

  cx =  @. focus * cos(ks)     + defocus * cosh(ks)
  cy =  @. focus * cosh(ks)    + defocus * cos(ks)
  sx =  @. focus * sincu(ks)   + defocus * sinhcu(ks)
  sy =  @. focus * sinhcu(ks)  + defocus * sincu(ks)
  sx2 = @. focus * sincu(2ks)  + defocus * sinhcu(2ks)
  sy2 = @. focus * sinhcu(2ks) + defocus * sincu(2ks)
  sxz = @. focus * sin(ks)^2   - defocus * sinh(ks)^2
  syz = @. focus * sinh(ks)^2  - defocus * sin(ks)^2

  @. vf.x  = vi.x  * cx + xp * s * sx
  @. vf.px = vi.px * cx - vi.x * p * k2 * s * sx
  @. vf.y  = vi.y  * cy + yp * s * sy
  @. vf.py = vi.py * cy + vi.y * p * k2 * s * sy
  @. vf.z  = (vi.z - (s / 4) * (xp^2 * (1 + sx2) + yp^2 * (1 + sy2)
                                + k2 * vi.x^2 * (1 - sx2) - k2 * vi.y^2 * (1 - sy2))
                   + (vi.x * xp * sxz - vi.y * yp * syz) / 2)
  @. vf.pz = vi.pz

  return vf
end # function trackQuadMx


"""
    trackQuadK!(vf, vi, s::Float64)

track "remaining part" of quadrupole, a position kick

### Implementation
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.
"""
function trackQuadK!(vf, vi, betgam_ref, s)
  tilde_m = 1 / betgam_ref  # mc^2 / p0·c
  beta_ref = sr_beta(betgam_ref)
  beta_ref = sr_beta(betgam_ref)
  gamsqr_ref = 1 + betgam_ref^2

  p    = @. 1 + vi.pz  # reduced momentum, P/P0 = 1 + δ
  ptr2 = @. vi.px^2 + vi.py^2
  ps   = @. sqrt(p^2 - ptr2)

  @. vf.x  = vi.x + s * vi.px / p * ptr2 / (ps * (p + ps))
  @. vf.y  = vi.y + s * vi.py / p * ptr2 / (ps * (p + ps))
  @. vf.z  = vi.z - s * ( (1.0 + vi.pz)
                         * (ptr2 - vi.pz * (2 + vi.pz) / gamsqr_ref)
                         / ( beta_ref * sqrt((1.0 + vi.pz)^2 + tilde_m^2) * ps
                             * (beta_ref * sqrt((1.0 + vi.pz)^2 + tilde_m^2) + ps)
                           ) - ptr2 / (2 * (1 + vi.pz)^2)
                        )
  @. vf.px = vi.px
  @. vf.py = vi.py
  @. vf.pz = vi.pz

  return vf
end # function trackQ!M::Quadrupole()


end # module MatrixKick

