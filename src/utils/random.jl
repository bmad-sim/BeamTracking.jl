"""
This function returns two Gaussian random numbers with 
mean 0 and standard deviations sigma1, sigma2.
"""
function gaussian_random(::CPU, sigma1, sigma2)
  return randn()*sigma1, randn()*sigma2
end


"""
See gaussian_random, but for SIMD vectors.
"""
function gaussian_random(::CPU, sigma1::SIMD.Vec, sigma2::SIMD.Vec)
  return SIMDMathFunctions.vmap((s1,s2)->(randn()*s1, randn()*s2), sigma1, sigma2)
end


"""
This function returns two Gaussian random numbers with 
mean 0 and standard deviations sigma1, sigma2 using a 
Box-Muller transform.

This was implemented because CUDA.randn has some horrible 
compiler bug, but CUDA.rand seems to be ok.
"""
function gaussian_random(::GPU, sigma1, sigma2)
  s, c = sincospi(2 * rand())
  t = sqrt(-2 * log(rand()))
  z0 = c*t*sigma1
  z1 = s*t*sigma2
  return z0, z1
end