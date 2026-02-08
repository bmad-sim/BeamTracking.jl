"""
Computes the mean coordinate vector and covariance matrix.
"""
function mean_and_cov(coords)
  v = coords.v 
  mu = SVector{6}(mean(v; dims=1))
  sigma = SMatrix{6,6}(cov(v; dims=1, corrected=false))
  return mu, sigma
end