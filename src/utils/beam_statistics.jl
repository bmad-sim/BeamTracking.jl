"""
Computes the mean coordinate vector and covariance matrix of v weighted by w.
"""
function mean_and_cov(v, w, ::CPU)
  W = sum(w)
  mean = SVector{6}((w' * v) ./ W)
  cov = SMatrix{6,6}((v' * (v .* w)) / W .- mean * mean')
  return mean, cov
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w, ::GPU)
  W = sum(w)
  mean = (w' * v) ./ W
  cov = SMatrix{6,6}((v' * (v .* w)) / W .- mean' * mean)
  return mean, cov
end


"""
Computes the mean coordinate vector and covariance matrix of v weighted uniformly.
"""
function mean_and_cov(v, w::Nothing, ::CPU)
  N = size(v, 1)
  mean = SVector{6}(sum(v; dims = 1) ./ N)
  cov = SMatrix{6,6}((v' * v) / N .- mean * mean')
  return mean, cov
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w::Nothing, ::GPU)
  N = size(v, 1)
  mean = sum(v; dims = 1) ./ N
  cov = SMatrix{6,6}((v' * v) / N .- mean' * mean)
  return mean, cov
end