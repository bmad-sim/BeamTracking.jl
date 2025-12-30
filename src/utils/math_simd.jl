# Straight from SIMD.jl:
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.ScalarTypes) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.ScalarTypes, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::T, v2::T) where {T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::T2, v2::SIMD.Vec{N, T}) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::T2) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
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
  #if isinf_real(x)
  #  return zero(x)
  #end

  threshold = 0.0004 # (120*eps(Float64))^(1/4)
  return vifelse(abs(x) > threshold, sin(x)/x, 1-x^2/6)
end


"""
    sinhcu(x)

Compute the unnormalized sinhc function ``\\operatorname{sinhcu}(x) = \\sinh(x) / (x)`` 
with accuracy near the origin.
"""
function sinhcu(x)
  #if isinf_imag(x)
  #  return zero(x)
  #end

  threshold = 0.0004 # (120*eps(Float64))^(1/4)
  return vifelse(abs(x) > threshold, sinh(x)/x, 1+x^2/6)
end


function atan2(y, x)
  return atan(y, x)
end


function atan2(y::SIMD.Vec{N, T}, x::SIMD.Vec{N, T}) where {N, T}
  arctan = atan(y/x)
  return vifelse(x > 0, arctan,
         vifelse((x < 0)  & (y >= 0),  arctan + SIMD.Vec{N, T}(T(pi)),
         vifelse((x < 0)  & (y < 0),   arctan - SIMD.Vec{N, T}(T(pi)),
         vifelse((x == 0) & (y > 0),   SIMD.Vec{N, T}(pi/2),
         vifelse((x == 0) & (y < 0),   SIMD.Vec{N, T}(-pi/2),
         vifelse((x == 0) & (y == 0),  SIMD.Vec{N, T}(0), 
         SIMD.Vec{N, T}(NaN)))))))
end


# Copy-pasted from sincc in bmad-ecosystem
function sincuc(x) 
  c0 = 1/6
  c1 = -1/120
  c2 = 1/5040
  c3 = -1/362880
  x2 = x^2
  return vifelse(abs(x) >= 0.1, (x-sin(x))/x^3, c0+x2*(c1+x2*(c2+x2*c3)))
end

"""
    one_cos_norm(x)

Function to calculate `(1 - cos(x)) / x^2` to machine precision.
This is usful if angle can be near zero where the direct evaluation of `(1 - cos(x))x^2` is inaccurate.
"""
one_cos_norm(x) = sincu(x/2)^2 / 2

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

"""
This function computes a Gaussian random number with mean zero and standard deviation sigma.
"""
function gaussian_random(sigma)
  return randn() * sigma 
end


"""
See gaussian_random, but for SIMD vectors.
"""
function gaussian_random(sigma::SIMD.Vec{N,T}) where {N,T}
  return vmap(gaussian_random, sigma)
end