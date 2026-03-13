"""
Computes the mean coordinate vector and covariance matrix.
"""
function mean_and_cov(v, ::CPU)
  N = size(v, 1)
  mu = SVector{6}(mean(v; dims=1))
  sigma = SMatrix{6,6}((v'*v)/N .- mu*mu')
  return mu, sigma
end

"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, ::GPU)
  N = size(v, 1)
  mu = SVector{6}(mean(v; dims=1))
  sigma = SMatrix{6,6}(Array((v'*v)/N .- mu*mu'))
  return mu, sigma
end