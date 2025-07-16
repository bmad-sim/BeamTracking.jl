#=

Utility functions and "fake" APC. These will be moved to 
AcceleratorSimUtils.jl in the end.

=#

#  Math =======================================================================
# u corresponds to unnormalized

# sinc/sincu is zero when the real part is Inf and imag is finite
isinf_real(x::Real) = isinf(x)
isinf_real(x::Number) = isinf(real(x)) && isfinite(imag(x))

# sinhc/sinhcu is zero when the imag part is Inf and real is finite
isinf_imag(x::Real) = false
isinf_imag(x::Number) = isfinite(real(x)) && isinf(imag(x))

# sincu copied from Boost library and correct limit behavior added
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinc.hpp
"""
    sincu(x)

Compute the unnormalized sinc function ``\\operatorname{sincu}(x) = \\sin(x) / (x)`` 
with accuracy near the origin.
"""
sincu(x) = _sinc(float(x))
function _sinc(x::Union{T,Complex{T}}) where {T}
    if isinf_real(x)
        return zero(x)
    end

    nrm = Base.Math.fastabs(x)
    if nrm >= 3.3*sqrt(sqrt(eps(T)))
        return sin(x)/x
    else
        # |x| < (eps*120)^(1/4)
        return 1 - x*x/6
    end
end

# sinhcu copied from Boost library and correct limit behavior added
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinhc.hpp

"""
    sinhcu(x)

Compute the unnormalized sinhc function ``\\operatorname{sinhcu}(x) = \\sinh(x) / (x)`` 
with accuracy accuracy near the origin.
"""
sinhcu(x) = _sinhcu(float(x))
function _sinhcu(x::Union{T,Complex{T}}) where {T}
    taylor_0_bound = eps(T)
    taylor_2_bound = sqrt(taylor_0_bound)
    taylor_n_bound = sqrt(taylor_2_bound)

    if isinf_imag(x) 
        return zero(x)
    end
    
    nrm = Base.Math.fastabs(x)

    if nrm >= taylor_n_bound || isnan(nrm)
        return sinh(x)/x
    else
        # approximation by taylor series in x at 0 up to order 0
        res = one(x)
        if nrm >= taylor_0_bound
            x2 = x*x
            # approximation by taylor series in x at 0 up to order 2
            res += x2/6
            if nrm >= taylor_2_bound
                # approximation by taylor series in x at 0 up to order 4
                res += (x2*x2)/120
            end
        end
        return res
    end
end

# Copied pasted from sincc in bmad-ecosystem
function sincuc(x) 
  if Base.Math.fastabs(x) < 0.1
    c0 = 1/6
    c1 = -1/120
    c2 = 1/5040
    c3 = -1/362880
    x2 = x^2
    y = c0 + x2 * (c1 + x2 * (c2 + x2 * c3))
  else
    y = (x - sin(x)) / x^3
  end
  return y
end

"""
returns instantaneous T_BMT quaternion in a straight coordinate system for particle `i` in a bunch
"""
function TBMT_quat(G, gx, gy, γ, Bρ, B::Vector{<:Number}, bv::Matrix{<:Number}, i)
  v = bv[i,:]
  rel_p = 1 + v[6]
  ps = sqrt(rel_p^2 - v[2]^2 - v[4]^2)
  β_hat = [v[2], v[4], ps] / rel_p
  B_para = (B' * β_hat) * β_hat
  B_perp = B - B_para
  # T-BMT precession in straight coordinate system
  Ω = -((1+G*γ) * B_perp .+ (1+G) * B_para) / (Bρ * ps)
  θ = norm(Ω)
  Ω /= (θ + (θ==0))  # if θ==0, denominator is 1, so Ω unchanged; else Ω /= θ
  return SVector(cos(θ/2), sin(θ/2)*Ω...)
end

function quat_mult!(q1, q2, i)
  # In-place quaternion multiplication: q2 := q1 * q2
  w1, x1, y1, z1 = q1
  w2, x2, y2, z2 = q2[i,:]
  q2[i,1] = w1*w2 - x1*x2 - y1*y2 - z1*z2
  q2[i,2] = w1*x2 + x1*w2 + y1*z2 - z1*y2
  q2[i,3] = w1*y2 - x1*z2 + y1*w2 + z1*x2
  q2[i,4] = w1*z2 + x1*y2 - y1*x2 + z1*w2
end

function quat_inv(q)
  w, x, y, z = q
  norm_sq = w^2 + x^2 + y^2 + z^2
  return ([w, -x, -y, -z] ./ norm_sq)
end

# Fake APC ====================================================================
const Q = 1.602176634e-19 # C
const C_LIGHT = 2.99792458e8 # m/s
const M_ELECTRON = 0.51099895069e6 # eV/c^2
const M_PROTON = 9.3827208943e8 # eV/c^2
const G_ELECTRON = 0.00115965218059
const G_PROTON = 1.79284734463

struct Species
  name::String
  mass::Float64   # in eV/c^2
  charge::Float64 # in Coulomb
  anomalous_magnetic_moment::Float64 # dimensionless, e.g. 0.001 for electron, 1.79285 for proton
end

const ELECTRON = Species("electron", M_ELECTRON,-1, G_ELECTRON)
const POSITRON = Species("positron", M_ELECTRON, 1, G_ELECTRON)

const PROTON = Species("proton", M_PROTON, 1, G_PROTON)
const ANTIPROTON = Species("antiproton", M_PROTON,-1, G_PROTON)


function Species(name)
  if name == "electron"
    return ELECTRON
  elseif name == "positron"
    return POSITRON
  elseif name == "proton"
    return PROTON
  elseif name == "ANTIPROTON"
    return ANTIPROTON
  else
    error("BeamTracking.jl's fake APC does not support species $name")
  end
end

massof(s::Species) = s.mass
chargeof(s::Species) = s.charge
anomalous_moment_of(s::Species) = s.anomalous_magnetic_moment

# Particle energy conversions =============================================================
calc_Brho(species::Species, E) = @FastGTPSA sqrt(E^2-massof(species)^2)/C_LIGHT/chargeof(species)
calc_E(species::Species, Brho) = @FastGTPSA sqrt(calc_p0c(species,Brho)^2 + massof(species)^2)
calc_gamma(species::Species, Brho) = @FastGTPSA sqrt((calc_p0c(species,Brho)/massof(species))^2+1)
calc_gamma(species::Species, Brho, δ) = @FastGTPSA sqrt((calc_p0c(species,Brho)*(1+δ))^2/massof(species)^2+1)

calc_p0c(species::Species, Brho) = @FastGTPSA Brho*C_LIGHT*chargeof(species)
calc_beta_gamma(species::Species, Brho) = @FastGTPSA calc_p0c(species,Brho)/massof(species)



#=


"""
    sr_gamma(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor ``\\gamma``.
"""
function sr_gamma(beta_gamma)
  return hypot(1, beta_gamma)
end



"""
    sr_gamma_m1(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor minus one, ``\\gamma - 1``.
"""
function sr_gamma_m1(beta_gamma)
  return beta_gamma^2 / (hypot(1, beta_gamma) + 1)
end


"""
    sr_beta(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the normalized velocity ``\\beta = v / c``.
"""
function sr_beta(beta_gamma)
  return beta_gamma / hypot(1, beta_gamma)
end


"""
    sr_pc(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the energy-equivalent momentum, ``pc``.
"""
function sr_pc(e_rest, beta_gamma)
  return e_rest * beta_gamma
end


"""
    sr_ekin(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the kinetic energy,
``E_\\text{kin} = mc^2(\\gamma - 1)``.
"""
function sr_ekin(e_rest, beta_gamma)
  return e_rest * sr_gamma_m1(beta_gamma)
end


"""
    sr_etot(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the total energy, ``E_\\tot = mc^2\\gamma``.
"""
function sr_etot(e_rest, beta_gamma)
  return e_rest * hypot(1, beta_gamma)
end


"""
    brho(e_rest, beta_gamma, ne = 1)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the magnetic rigidity, ``B\\rho = p / q``.

DTA: Need to handle energy units other than ``\\mathrm{eV}``..
"""
function brho(e_rest, beta_gamma, ne = 1)
  return (sr_pc(e_rest, beta_gamma) / (ne * C_LIGHT))
end
## If given ``E_\text{kin}`` instead of ``\beta\gamma``,
## use the following:
#
#function sr_gamma(e_rest, e_kin)
#  return e_kin / e_rest + 1
#end
#
#function sr_gamma_m1(e_rest, e_kin)
#  return e_kin / e_rest
#end
#
#function sr_beta_gamma(e_rest, e_kin)
#  return sqrt(e_kin / e_rest * (e_kin / e_rest + 2))
#end
#
#function sr_beta(e_rest, e_kin)
#  return sr_beta_gamma(e_rest, e_kin) / sr_gamma(e_rest, e_kin)
#end
#
#function sr_pc(e_rest, e_kin)
#  #return sqrt(e_kin * (e_kin + 2e_rest))
#  return e_rest * sr_beta_gamma(e_rest, e_kin)
#end
#
#function sr_etot(e_rest, e_kin)
#  return e_rest + e_kin
#end
#
#function brho(e_rest, e_kin, ne = 1)
#  return sr_pc(e_rest, e_kin) / (ne * clight)
#end

"""
    sincu(z)

## Description
Compute the unnormalized sinc function, ``\\operatorname{sincu}(z) = \\sin(z) / z``,
with a correct treatment of the removable singularity at the origin.

### Implementation
The function ``\\sin(z) / z = 1 - z^2 / 3! + z^4 / 5! - z^6 / 7! + \\cdots``.
For real values of ``z \\in (-1,1)``, one can truncate this series just before
the ``k^\\text{th}`` term, ``z^{2k} / (2k+1)!``, and the alternating nature of
the series guarantees the error does not exceed the value of this first truncated term.
It follows that if ``z^2 / 6 < \\varepsilon_\\text{machine}``, simply truncating
the series to 1 yields machine precision accuracy near the origin.
And outside that neighborhood, one simply computes ``\\sin(z) / z``.
On the otherhand, if we allow for complex values, one can no longer assume the
series alternates, and one must use a more conservative threshold.
Numerical exploration suggests that for ``|z| < 1`` the sum of terms starting
at the ``k^\\text{th}`` term is bounded by ``|z|^{2k} / (2k)!``.
It then follows that one should use ``|z|^2 / 2 < \\varepsilon_\\text{machine}``
as the criterion for deciding to truncate the series to 1 near the origin.
"""
function sincu(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sin(z) / z
  end
end

"""
    sinhcu(z)

## Description
Compute the hyperbolic sinc function, ``\\operatorname{sinhcu}(z) = \\operatorname{sinh}(z) / z``,
with a correct treatment of the removable singularity at the origin.

### Implementation
See sincu for notes about implementation.
"""
function sinhcu(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sinh(z) / z
  end
end

"""
    get_work(bunch::Bunch, ::Val{N}) where N -> work

Returns a tuple of `N` arrays of type `eltype(Bunch.v.x)` and 
length `length(Bunch.v.x)` which may be used as temporaries.

### Arguments
- `bunch`     -- Bunch to extract types and number of particles from
- `::Val{N}` -- Number of `N` temporary arrays desired
"""
function get_work(bunch::Bunch, ::Val{N}) where {N}
  sample = first(bunch.v.x)
  T = typeof(sample)
  N_particle = length(bunch.v.x)

  # Instead of using zeros, we do this to ensure 
  # same GTPSA descriptor if T isa TPS.
  return ntuple(Val{N}()) do t
    r = Vector{T}(undef, N_particle)
    for idx in eachindex(r)
      r[idx] = zero(sample)
    end
    r
  end
end

=#