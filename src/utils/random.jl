"""
This function returns two Gaussian random numbers with 
mean 0 and standard deviations sigma1, sigma2.
"""
function gaussian_random(::CPU, sigma1::T1, sigma2::T2) where {T1,T2}
  return randn(T1)*sigma1, randn(T2)*sigma2
end


"""
See gaussian_random, but for SIMD vectors.
"""
function gaussian_random(::CPU, sigma1::SIMD.Vec{N,T1}, sigma2::SIMD.Vec{N,T2}) where {N,T1,T2}
  return SIMDMathFunctions.vmap((s1,s2)->(randn(T1)*s1, randn(T2)*s2), sigma1, sigma2)
end


"""
This function returns two Gaussian random numbers with 
mean 0 and standard deviations sigma1, sigma2 using a 
Box-Muller transform.

This was implemented because CUDA.randn has some horrible 
compiler bug, but CUDA.rand seems to be ok.
"""
function gaussian_random(::GPU, sigma1::T1, sigma2::T2) where {T1, T2}
  s, c = sincospi(2 * rand(T1))
  t = sqrt(-2 * log(rand(T2)))
  z0 = c*t*sigma1
  z1 = s*t*sigma2
  return z0, z1
end