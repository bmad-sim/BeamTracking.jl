"""
Computes the mean coordinate vector and covariance matrix of v weighted uniformly.
"""
function mean_and_cov(v, w::Union{Number,Nothing}, ::CPU)
  N = size(v, 1)
  mu = mean(v; dims = 1)
  v_c = v .- mu
  sigma = (v_c' * v_c) ./ N
  return SVector{6}(mu), SMatrix{6,6}(sigma)
end


"""
Computes the mean coordinate vector and covariance matrix of v weighted by w.
"""
function mean_and_cov(v, w, ::CPU)
  W = sum(w)
  mu = (w' * v) ./ W
  v_c = v .- mu
  v_c = v_c .* sqrt.(w)
  sigma = (v_c' * v_c) ./ W
  return SVector{6}(mu), SMatrix{6,6}(sigma)
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w::Union{Number,Nothing}, ::GPU)
  N = size(v, 1)
  mu = mean(v; dims = 1)
  v_c = v .- mu
  sigma = Array((v_c' * v_c) ./ N)
  return SVector{6}(Array(mu)), SMatrix{6,6}(sigma)
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w, ::GPU)
  W = sum(w)
  mu = (w' * v) ./ W
  v_c = v .- mu
  v_c = v_c .* sqrt.(w)
  sigma = Array((v_c' * v_c) ./ W)
  return SVector{6}(Array(mu)), SMatrix{6,6}(sigma)
end
