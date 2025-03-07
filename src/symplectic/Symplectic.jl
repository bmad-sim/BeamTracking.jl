module Symplectic
using ..BeamTracking: Coords, Beam
using ..GTPSA: @FastGTPSA!

export track!


struct Drift{T}
  L::T  # drift length / m
end

struct Quadrupole{T}
  L::T   # quadrupole length / m
  B1::T  # quadrupole gradient / (T·m^-1)
end


"""
    track!(ele::Symplectic.Drift, beamf::Beam, beami::Beam) -> beamf

## Description
This function performs exact, hence symplectic, tracking of a `Beam` through a drift.

### Arguments
  - `ele`   -- element of type `Symplectic.Drift`
  - `beamf` -- final `Beam`
  - `beami` -- initial `Beam`
"""
function track!(ele::Drift, beamf::Beam, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  @assert !(beamf.species === beami.species) "Input species must be equal to output!"
  L = ele.L
  zi = beami.z
  zf = beamf.z

  tilde_m = 1 / beami.beta_gamma_0
  beta_0 = sr_beta(beami.beta_gamma_0)

  @. ps = sqrt((1.0 + zi.pz)^2 - (zi.px^2 + zi.py^2))
  @. et = sqrt((1.0 + zi.pz)^2 + tilde_m^2)

  @. zf.x = zi.x + zi.px * L / ps
  @. zf.px = zi.px
  @. zf.y = zi.y + zi.py * L / ps
  @. zf.py = zi.py
  @. zf.z = zi.z - (1.0 + zi.pz) * (L / ps - L / (beta_0 * et))
  @. zf.pz = zi.pz

  return beamf
end # function track!(::Drift, ::Beam, ::Beam)


"""
    track!(ele::Quadrupole, beamf::Beam, beami::Beam)

## Description
track quadrupole


### Arguments`
  - `ele`   -- element of type `Symplectic.Quadrupole`
  - `beamf` -- final `Beam`
  - `beami` -- initial `Beam`

### Implementation
This integrator uses the so-called Matrix-Kick-Matrix method to implement
an integrator accurate though second-order in the integration step-size.
"""
function track!(ele::Symplectic.Quadrupole, beamf::Beam, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  L = ele.L

  # κ^2 (kappa-squared) = (q g / P0) / (1 + δ)
  # numerator of κ^2
  k2_num = chargeof(species) * ele.B1 / ele.P0

  trackQuadMx!(beamf, beami, k2_num, L / 2)
  trackQuadK!( beami, beamf, L)
  trackQuadMx!(beamf, beami, k2_num, L / 2)

  return beamf
end # function track!(::Quadrupole)


"""
    trackQuadMx!(beamf::Beam, beami::Beam, k2_num::Float64, s::Float64)

track "matrix part" of quadrupole
"""
function trackQuadMx!(beamf::Beam, beami::Beam, k2_num, s)
  zi = beami.z
  zf = beamf.z
  focus   = k2_num >= 0  # horizontally focusing
  defocus = k2_num <  0  # horizontally defocusing

  @. p = 1 + zi.pz    # reduced momentum, P/P0 = 1 + δ
  @. k2 = k2_num / p  # κ^2 for each particle
  @. ks = sqrt(abs(k2)) * s  # |κ|s
  @. xp = px / p  # x'
  @. yp = py / p  # y'

  @. cx =  focus * cos(ks)    + defocus * cosh(ks)
  @. cy =  focus * cosh(ks)   + defocus * cos(ks)
  @. sx =  focus * sincu(ks)  + defocus * sinch(ks)
  @. sy =  focus * sinch(ks)  + defocus * sincu(ks)
  @. sx2 = focus * sincu(2ks) + defocus * sinch(2ks)
  @. sy2 = focus * sinch(2ks) + defocus * sincu(2ks)
  @. sxz = focus * sin(ks)^2  + defocus * sinh(ks)^2
  @. syz = focus * sinh(ks)^2 + defocus * sin(ks)^2

  @. zf.x  = zi.x  * cx + xp * s * sx
  @. zf.px = zi.px * cx - k2 * p * zi.x * s * sx
  @. zf.y  = zi.y  * cy + yp * s * sy
  @. zf.py = zi.py * cy + k2 * p * zi.y * s * sy
  @. zf.z  = (zi.z - (s / 4) * (xp^2 * (1 + sx2) + yp^2 * (1 + sy2)
                                + k2 * zi.x^2 * (1 - sx2) - k2 * zi.y^2 * (1 - sy2))
                   + (zi.x * xp * sxz - zi.y * yp * syz) / 2.)
  @. zf.pz = zi.pz

  return beamf
end # function trackQuadMx


"""
    trackQuadK!(beamf::Beam, beami::Beam, s::Float64)

track "remaining part" of quadrupole, a position kick

### Implementation
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.
"""
function trackQuadK!(beamf::Beam, beami::Beam, s::Float64)
  zi = beami.z
  zf = beamf.z

  tilde_m = 1 / beami.beta_gamma_0
  beta_0 = sr_beta(beami.beta_gamma_0)

  @. p = 1 + zi.pz  # reduced momentum, P/P0 = 1 + δ
  @. ptr2 = px^r + py^2
  @. ps = sqrt(p^2 - ptr2)

  @. zf.x  = zi.x + s * zi.px / p * ptr2 / (ps (p + ps))
  @. zf.y  = zi.y + s * zi.py / p * ptr2 / (ps (p + ps))
  @. zf.z  = zi.z - s * (p / ps - (1/2) * ptr2 / p^2
                         - p / (β0 * sqrt(p^2 + tilde_m^2)))
  @. zf.px = zi.px
  @. zf.py = zi.py
  @. zf.pz = zi.pz

  return beamf
end # function trackQ!M::Quadrupole()


end # module Symplectic

