"""
Computes the mean coordinate vector and covariance matrix of the rows of v selected by the
mask `alive`, weighted uniformly.

Unselected rows are excluded with `ifelse` rather than multiplied by zero, which would
propagate NaN or Inf coordinates (0*NaN is NaN).
"""
function mean_and_cov(v, w::Union{Number,Nothing}, alive, ::CPU)
  n = count(alive)
  mu = reshape([sum(i -> @inbounds(ifelse(alive[i], v[i,j], zero(eltype(v)))), axes(v, 1)) for j in axes(v, 2)] ./ n, 1, :)
  v_c = ifelse.(alive, v .- mu, zero(eltype(mu)))
  sigma = (v_c' * v_c) ./ n
  return SVector{6}(mu), SMatrix{6,6}(sigma)
end


"""
Computes the mean coordinate vector and covariance matrix of the rows of v selected by the
mask `alive`, weighted by w.
"""
function mean_and_cov(v, w, alive, ::CPU)
  W = sum(i -> @inbounds(ifelse(alive[i], w[i], zero(eltype(w)))), axes(v, 1))
  mu = reshape([sum(i -> @inbounds(ifelse(alive[i], w[i]*v[i,j], zero(eltype(w))*zero(eltype(v)))), axes(v, 1)) for j in axes(v, 2)] ./ W, 1, :)
  # Unit weights at the unselected rows keep sqrt away from zero, which GTPSA rejects
  T = promote_type(eltype(mu), typeof(sqrt(one(eltype(w)))))
  v_c = ifelse.(alive, (v .- mu) .* sqrt.(ifelse.(alive, w, one(eltype(w)))), zero(T))
  sigma = (v_c' * v_c) ./ W
  return SVector{6}(mu), SMatrix{6,6}(sigma)
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w::Union{Number,Nothing}, alive, ::GPU)
  n = count(alive)
  mu = sum(ifelse.(alive, v, zero(eltype(v))); dims = 1) ./ n
  v_c = ifelse.(alive, v .- mu, zero(eltype(mu)))
  sigma = Array((v_c' * v_c) ./ n)
  return SVector{6}(Array(mu)), SMatrix{6,6}(sigma)
end


"""
See mean_and_cov, but on the GPU.
"""
function mean_and_cov(v, w, alive, ::GPU)
  w_alive = ifelse.(alive, w, zero(eltype(w)))
  W = sum(w_alive)
  mu = (w_alive' * ifelse.(alive, v, zero(eltype(v)))) ./ W
  T = promote_type(eltype(mu), typeof(sqrt(one(eltype(w)))))
  v_c = ifelse.(alive, (v .- mu) .* sqrt.(ifelse.(alive, w, one(eltype(w)))), zero(T))
  sigma = Array((v_c' * v_c) ./ W)
  return SVector{6}(Array(mu)), SMatrix{6,6}(sigma)
end
