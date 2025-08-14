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
    sincus(x)

Compute the unnormalized sinc square-root function 
``\\operatorname{sincus}(x) = \\sin(\\sqrt(x)) / (\\sqrt(x))`` 
with accuracy near the origin.
"""
sincus(x) = _sincus(float(x))
function _sincus(x::Union{T,Complex{T}}) where {T}
    nrm = Base.Math.fastabs(x)
    if nrm >= 109.018*eps(T)^(1/11)
        return sin(sqrt(x))/(sqrt(x))
    else
        c0 = 1
        c1 = -1/6
        c2 = 1/120
        c3 = -1/5040
        c4 = 1/362880
        c5 = -1/39916800
        c6 = 1/6227020800
        c7 = -1/1307674368000
        c8 = 1/355687428096000
        c9 = -1/12164510040883200000
        c10 = 1/5109094217170944000000
        return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*
        (c9+x*c10))))))))))
    end
end

"""
    coss(x)

Compute the cos square-root function 
``\\operatorname{coss}(x) = \\cos(\\sqrt(x))`` 
with differentiability near the origin.
"""
coss(x) = _coss(float(x))
function _coss(x::Union{T,Complex{T}}) where {T}
    nrm = Base.Math.fastabs(x)
    if nrm >= 81.9796*eps(T)^(1/11)
        return cos(sqrt(x))
    else
        c0 = 1
        c1 = -1/2
        c2 = 1/24
        c3 = -1/720
        c4 = 1/40320
        c5 = -1/3628800
        c6 = 1/479001600
        c7 = -1/87178291200
        c8 = 1/20922789888000
        c9 = -1/6402373705728000
        c10 = 1/2432902008176640000
        return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*
        (c9+x*c10))))))))))
    end
end

@inline function expq(v)
  """
  This function computes exp(i v⋅σ) as a quaternion, where σ is the 
  vector of Pauli matrices.
  """
  @FastGTPSA begin
    n2 = v[1]^2 + v[2]^2 + v[3]^2
    c = coss(n2)
    s = sincus(n2)
  v2 = s * v
  end
  return SA[-c, v2[1], v2[2], v2[3]]
end

@inline function quat_mul(q1, q2)
  """
  Returns q1 * q2.
  """
  @FastGTPSA begin
    a1, b1, c1, d1 = q1[Q0], q1[QX], q1[QY], q1[QZ]
    a2, b2, c2, d2 = q2[Q0], q2[QX], q2[QY], q2[QZ]
    a3 = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b3 = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c3 = a1*c2 - b1*d2 + c1*a2 + d1*b2
    d3 = a1*d2 + b1*c2 - c1*b2 + d1*a2
  end
  return SA[a3 b3 c3 d3]
end

# Particle energy conversions =============================================================
R_to_E(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT*chargeof(species))^2 + massof(species)^2)
E_to_R(species::Species, E) = @FastGTPSA massof(species)*sinh(acosh(E/massof(species)))/C_LIGHT/chargeof(species) 
pc_to_R(species::Species, pc) = @FastGTPSA pc/C_LIGHT/chargeof(species)

R_to_gamma(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT/massof(species))^2+1)
R_to_pc(species::Species, R) = @FastGTPSA R*chargeof(species)*C_LIGHT
R_to_beta_gamma(species::Species, R_ref) = @FastGTPSA R_ref*chargeof(species)*C_LIGHT/massof(species)

# Fake APC because APC is not working for now :(
function anom(species::Species) 
  if nameof(species) == "electron"
    return 0.00115965218046
  elseif nameof(species) == "positron"
    return 0.0011596521735304233
  elseif nameof(species) == "proton"
    return 1.7928473446300592
  elseif nameof(species) == "proton"
    return 1.7928473446300592
  else
    error("Your species is not in fake APC yet")
  end
end

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
    R_ref(e_rest, beta_gamma, ne = 1)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the magnetic R_ref, ``B\\rho = p / q``.

DTA: Need to handle energy units other than ``\\mathrm{eV}``..
"""
function R_ref(e_rest, beta_gamma, ne = 1)
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
#function R_ref(e_rest, e_kin, ne = 1)
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
