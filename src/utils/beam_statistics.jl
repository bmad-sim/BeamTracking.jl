"""
Computes the mean vector and covariance matrix.
"""
function mean_and_cov(coords)
  v = coords.v 
  mu = mean(v; dims=1)
  sigma = cov(v; dims=1, corrected=false)
  return mu, sigma
end