# Straight from SIMD.jl:
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.ScalarTypes) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.ScalarTypes, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::T, v2::T) where {T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::T2, v2::SIMD.Vec{N, T}) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::T2) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
# One that probably should be in SIMD but isn't (?)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::T, v2::T) where {N,T<:SIMD.ScalarTypes} = vifelse(v, SIMD.Vec{N,T}(v1), SIMD.Vec{N,T}(v2))
# Fallback for type unstable:
@inline vifelse(v::Union{Bool,SIMD.Vec{N, Bool}}, v1, v2) where {N} = ifelse(v, v1, v2)


# sinc/sincu is zero when the real part is Inf and imag is finite
isinf_real(x::Real) = isinf(x)
isinf_real(x::Number) = isinf(real(x)) && isfinite(imag(x))

# sinhc/sinhcu is zero when the imag part is Inf and real is finite
isinf_imag(x::Real) = false
isinf_imag(x::Number) = isfinite(real(x)) && isinf(imag(x))


"""
    sincu(x)

Compute the unnormalized sinc function ``\\operatorname{sincu}(x) = \\sin(x) / (x)`` 
with accuracy near the origin.
"""
function sincu(x)
  # Taylor series through x^8 for |x| < 0.1: its error x^10/39916800 is below eps(Float64)/2,
  # and its derivatives (ForwardDiff) stay accurate where those of sin(x)/x cancel
  x2 = x*x
  s = 1 - x2/6*(1 - x2/20*(1 - x2/42*(1 - x2/72)))
  return vifelse(abs(x) < one(x)/10, s, sin(x)/x)
end


"""
    sinhcu(x)

Compute the unnormalized sinhc function ``\\operatorname{sinhcu}(x) = \\sinh(x) / (x)`` 
with accuracy near the origin.
"""
function sinhcu(x)
  # Taylor series through x^8 for |x| < 0.1, as in sincu
  x2 = x*x
  s = 1 + x2/6*(1 + x2/20*(1 + x2/42*(1 + x2/72)))
  return vifelse(abs(x) < one(x)/10, s, sinh(x)/x)
end


function atan2(y, x)
  return atan(y, x)
end


@generated function atan2(y::SIMD.Vec{N, T}, x::SIMD.Vec{N, T}) where {N, T}
  p = T(pi)
  hp = T(pi/2)
  n = T(NaN)
  return quote
    arctan = atan(y/x)
    return vifelse(x > 0, arctan,
           vifelse((x < 0)  & (y >= 0),  arctan + SIMD.Vec{N, T}($p),
           vifelse((x < 0)  & (y < 0),   arctan - SIMD.Vec{N, T}($p),
           vifelse((x == 0) & (y > 0),   SIMD.Vec{N, T}($hp),
           vifelse((x == 0) & (y < 0),   SIMD.Vec{N, T}(-$hp),
           vifelse((x == 0) & (y == 0),  SIMD.Vec{N, T}(0), 
           SIMD.Vec{N, T}($n)))))))
  end
end


# The exact bend evaluates asin, sincos, atan2 and sin per particle step, which do not vectorize
# (SIMD Vecs run them lane by lane). When the argument of every particle is small, a Taylor series
# is used instead (truncation below 1e-16 relative; one compare, plus one reduction for a Vec).
# Only for floats and Vecs of them: TPS and Duals keep the functions, whose derivatives of every
# order are exact (a series truncated at x^11 gives wrong high-order TPS coefficients)
const SeriesTypes = Union{AbstractFloat, SIMD.Vec{<:Any,<:AbstractFloat}}

asin_s(x) = asin(x)
function asin_s(x::SeriesTypes)
  # Through x^11 for |x| < 0.05
  all(abs(x) < one(x)/20) || return asin(x)
  x2 = x*x
  return x*(1 + x2*(one(x)/6 + x2*(3*one(x)/40 + x2*(5*one(x)/112 + x2*(35*one(x)/1152 +
         x2*(63*one(x)/2816))))))
end

sincos_s(x) = sincos(x)
function sincos_s(x::SeriesTypes)
  # Through x^9 (sin) and x^10 (cos) for |x| < 0.1
  all(abs(x) < one(x)/10) || return sincos(x)
  x2 = x*x
  s = x*(1 - x2*(one(x)/6)*(1 - x2*(one(x)/20)*(1 - x2*(one(x)/42)*(1 - x2*(one(x)/72)))))
  c = 1 - x2*(one(x)/2)*(1 - x2*(one(x)/12)*(1 - x2*(one(x)/30)*(1 - x2*(one(x)/56)*
      (1 - x2*(one(x)/90)))))
  return s, c
end

sin_s(x) = sin(x)
sin_s(x::SeriesTypes) = all(abs(x) < one(x)/10) ? first(sincos_s(x)) : sin(x)

atan2_s(y, x) = atan2(y, x)
function atan2_s(y::T, x::T) where {T<:SeriesTypes}
  # atan(y/x) through (y/x)^11 for x > 0 and |y/x| < 0.05
  r = y/x
  all((x > 0) & (abs(r) < one(r)/20)) || return atan2(y, x)
  r2 = r*r
  return r*(1 - r2*(one(r)/3 - r2*(one(r)/5 - r2*(one(r)/7 - r2*(one(r)/9 - r2*(one(r)/11))))))
end


# Copy-pasted from sincc in bmad-ecosystem
function sincuc(x) 
  c0 = 1//6
  c1 = -1//120
  c2 = 1//5040
  c3 = -1//362880
  x2 = x^2
  return vifelse(abs(x) >= 1//10, (x-sin(x))/x^3, c0+x2*(c1+x2*(c2+x2*c3)))
end

"""
    one_cos_norm(x)

Function to calculate `(1 - cos(x)) / x^2` to machine precision.
This is useful if angle can be near zero where the direct evaluation of 
`(1 - cos(x)) / x^2` is inaccurate.
"""
one_cos_norm(x) = sincu(x/2)^2/2

#=
"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
@inline function bessel01_RF(x)
  threshold = 2.9e-7 # sqrt(64*eps(Float64))
  sq = sqrt(x)
  b0_out = besselj(0, sq)
  b1 = besselj(1, sq)
  b1 = b1/sq
  b1_out = ifelse(x > threshold, b1, 1/2-x/16)
  return b0_out, b1_out
end


"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
@inline function bessel01_RF(x::SIMD.Vec{N, T}) where {N, T}
  threshold = 2.9e-7 # sqrt(64*eps(Float64))
  sq = sqrt(x)
  b0_out = SIMDMathFunctions.vmap(besselj0, sq)
  b1 = SIMDMathFunctions.vmap(besselj1, sq)
  b1 = b1/sq
  b1_out = vifelse(x > threshold, b1, 1/2-x/16)
  return b0_out, b1_out
end


"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
function bessel01_RF(x::TPS{T}) where {T}
  ε = eps(T)
  N_max = 100
  N = 1
  conv0 = false
  conv1 = false
  y = one(x)
  prev0 = one(x)
  prev1 = one(x)
  result0 = one(x)
  result1 = one(x)/2
  @FastGTPSA begin
    while !(conv0 && conv1) && N < N_max
      y = -y*x/(4*N*N)
      result0 = result0 + y
      result1 = result1 + y/(2*(N + 1))
      N += 1
      if normTPS(result0 - prev0) < ε
        conv0 = true
      end
      if normTPS(result1 - prev1) < ε
        conv1 = true
      end
      prev0 = result0
      prev1 = result1
    end
  end
  if N == N_max
    @warn "bessel01_RF convergence not reached in $N_max iterations"
  end
  return result0, result1
end
=#