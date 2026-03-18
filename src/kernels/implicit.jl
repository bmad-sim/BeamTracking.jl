@inline norm_sq(v::SVector{3,T}) where {T} = muladd(v[1], v[1], muladd(v[2], v[2], v[3] * v[3]))

@inline function ∇q_H(A::F, Jac::J, q::SVector{3,T}, p::SVector{3,T}, s, g, beta_0, tilde_m2) where {F,J,T}
  p_x, p_y, p_z = p[1], p[2], p[3]
  x, y, z = q[1], q[2], q[3]
  rel_p = 1 + p_z
  tilde_E2 = muladd(rel_p, rel_p, tilde_m2)

  beta_inv = beta_0 * sqrt(tilde_E2)
  beta = rel_p / beta_inv
  t = z / (beta * C_LIGHT)

  _A = A(x, y, s, t)
  _J = Jac(x, y, s, t)

  kappa = muladd(g, x, 1)
  p_kin_x = p_x - _A[2]
  p_kin_y = p_y - _A[3]
  p_kin_2 = muladd(p_kin_x, p_kin_x, p_kin_y * p_kin_y)
  p_s2 = muladd(rel_p, rel_p, -p_kin_2)
  p_s = sqrt(vifelse(p_s2 > 0, p_s2, one(p_s2)))
  scale = kappa / p_s

  dH_dqx = -_J[4, 1] - _J[1, 1] - scale * muladd(p_kin_x, _J[2, 1], p_kin_y * _J[3, 1]) - g * p_s
  dH_dqy = -_J[4, 2] - _J[1, 2] - scale * muladd(p_kin_x, _J[2, 2], p_kin_y * _J[3, 2])
  dH_dqz = -_J[4, 3] - _J[1, 3] - scale * muladd(p_kin_x, _J[2, 3], p_kin_y * _J[3, 3])

  return SVector{3,T}(dH_dqx, dH_dqy, dH_dqz)
end

@inline function ∇p_H(A::F, q::SVector{3,T}, p::SVector{3,T}, s, g, beta_0, tilde_m2) where {F,T}
  p_x, p_y, p_z = p[1], p[2], p[3]
  x, y, z = q[1], q[2], q[3]
  rel_p = 1 + p_z

  tilde_E2 = muladd(rel_p, rel_p, tilde_m2)
  beta_inv = beta_0 * sqrt(tilde_E2)
  beta = rel_p / beta_inv
  t = z / (beta * C_LIGHT)

  A_val = A(x, y, s, t)
  kappa = muladd(g, x, 1)

  p_kin_x = p_x - A_val[2]
  p_kin_y = p_y - A_val[3]
  p_kin_2 = muladd(p_kin_x, p_kin_x, p_kin_y * p_kin_y)
  p_s2 = muladd(rel_p, rel_p, -p_kin_2)
  p_s = sqrt(vifelse(p_s2 > 0, p_s2, one(p_s2)))
  scale = kappa / p_s

  dH_dpx = scale * p_kin_x
  dH_dpy = scale * p_kin_y
  dH_dpz = beta - scale * rel_p

  return SVector{3,T}(dH_dpx, dH_dpy, dH_dpz)
end

@inline function solve_3x3(A::SMatrix{3,3,T}, b::SVector{3,T}) where {T}
  c11 = A[2, 2] * A[3, 3] - A[2, 3] * A[3, 2]
  c12 = A[2, 3] * A[3, 1] - A[2, 1] * A[3, 3]
  c13 = A[2, 1] * A[3, 2] - A[2, 2] * A[3, 1]
  c21 = A[1, 3] * A[3, 2] - A[1, 2] * A[3, 3]
  c22 = A[1, 1] * A[3, 3] - A[1, 3] * A[3, 1]
  c23 = A[1, 2] * A[3, 1] - A[1, 1] * A[3, 2]
  c31 = A[1, 2] * A[2, 3] - A[1, 3] * A[2, 2]
  c32 = A[1, 3] * A[2, 1] - A[1, 1] * A[2, 3]
  c33 = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
  inv_det = one(T) / muladd(A[1, 1], c11, muladd(A[1, 2], c12, A[1, 3] * c13))
  return SVector{3,T}(
    muladd(c11, b[1], muladd(c21, b[2], c31 * b[3])) * inv_det,
    muladd(c12, b[1], muladd(c22, b[2], c32 * b[3])) * inv_det,
    muladd(c13, b[1], muladd(c23, b[2], c33 * b[3])) * inv_det
  )
end




@inline function hessian_pq(q::SVector{3,T}, p, A, Jac, s, g, beta_0, tilde_m2) where {T}
  x, y, z = q
  p_x, p_y, p_z = p
  rel_p = 1 + p_z
  tilde_E2 = muladd(rel_p, rel_p, tilde_m2)
  beta_inv = beta_0 * sqrt(tilde_E2)
  beta = rel_p / beta_inv
  t = z / (beta * C_LIGHT)

  _A = A(x, y, s, t)
  _J = Jac(x, y, s, t)
  p_kin_x = p_x - _A[2]
  p_kin_y = p_y - _A[3]
  p_s2 = muladd(rel_p, rel_p, -(muladd(p_kin_x, p_kin_x, p_kin_y * p_kin_y)))
  p_s = sqrt(vifelse(p_s2 > 0, p_s2, one(p_s2)))
  inv_ps = one(T) / p_s
  inv_ps2 = inv_ps * inv_ps
  inv_ps3 = inv_ps2 * inv_ps
  kappa = muladd(g, x, one(T))
  dkappa_dx = g

  function col(j)
    dpxk = -_J[2, j]
    dpyk = -_J[3, j]
    dps = (p_kin_x * _J[2, j] + p_kin_y * _J[3, j]) * inv_ps
    dkappa = (j == 1) ? dkappa_dx : zero(T)
    dscale = dkappa * inv_ps - kappa * dps * inv_ps2
    scale = kappa * inv_ps
    ddpx = dscale * p_kin_x + scale * dpxk
    ddpy = dscale * p_kin_y + scale * dpyk
    dddp = -dscale * rel_p
    return SVector{3,T}(ddpx, ddpy, dddp)
  end

  c1 = col(1)
  c2 = col(2)
  c3 = col(3)
  return SMatrix{3,3,T}(c1[1], c1[2], c1[3], c2[1], c2[2], c2[3], c3[1], c3[2], c3[3])
end

@inline function hessian_qp(q, p::SVector{3,T}, A, Jac, s, g, beta_0, tilde_m2) where {T}
  x, y, z = q
  p_x, p_y, p_z = p
  rel_p = 1 + p_z
  tilde_E2 = muladd(rel_p, rel_p, tilde_m2)
  beta_inv = beta_0 * sqrt(tilde_E2)
  beta = rel_p / beta_inv
  t = z / (beta * C_LIGHT)

  _A = A(x, y, s, t)
  _J = Jac(x, y, s, t)
  p_kin_x = p_x - _A[2]
  p_kin_y = p_y - _A[3]
  p_s2 = muladd(rel_p, rel_p, -(muladd(p_kin_x, p_kin_x, p_kin_y * p_kin_y)))
  p_s = sqrt(vifelse(p_s2 > 0, p_s2, one(p_s2)))
  inv_ps = one(T) / p_s
  inv_ps2 = inv_ps * inv_ps
  inv_ps3 = inv_ps2 * inv_ps
  kappa = muladd(g, x, one(T))
  scale = kappa * inv_ps

  dps_dpx = -p_kin_x * inv_ps
  dps_dpy = -p_kin_y * inv_ps
  dps_dpz = rel_p * inv_ps
  dscale_dpx = -kappa * dps_dpx * inv_ps2
  dscale_dpy = -kappa * dps_dpy * inv_ps2
  dscale_dpz = -kappa * dps_dpz * inv_ps2

  function col(dpidx)
    dscale = dpidx == 1 ? dscale_dpx : (dpidx == 2 ? dscale_dpy : dscale_dpz)
    dpxk = dpidx == 1 ? one(T) : zero(T)
    dpyk = dpidx == 2 ? one(T) : zero(T)
    dps = dpidx == 1 ? dps_dpx : (dpidx == 2 ? dps_dpy : dps_dpz)
    dT1 = muladd(dpxk, _J[2, 1], dpyk * _J[3, 1])
    dT2 = muladd(dpxk, _J[2, 2], dpyk * _J[3, 2])
    dT3 = muladd(dpxk, _J[2, 3], dpyk * _J[3, 3])
    T1 = muladd(p_kin_x, _J[2, 1], p_kin_y * _J[3, 1])
    T2 = muladd(p_kin_x, _J[2, 2], p_kin_y * _J[3, 2])
    T3 = muladd(p_kin_x, _J[2, 3], p_kin_y * _J[3, 3])
    ddq1 = -(dscale * T1 + scale * dT1) - g * dps
    ddq2 = -(dscale * T2 + scale * dT2)
    ddq3 = -(dscale * T3 + scale * dT3)
    return SVector{3,T}(ddq1, ddq2, ddq3)
  end

  c1 = col(1)
  c2 = col(2)
  c3 = col(3)
  return SMatrix{3,3,T}(c1[1], c1[2], c1[3], c2[1], c2[2], c2[3], c3[1], c3[2], c3[3])
end

#= =============================================================
# Can be done with Descriptor(6,1) since Hamiltonian is known
# Don't want to allocate new descriptor for each step
#
# @inline function hessian_pq_source(q::SVector{3,T}, p::SVector{3,T}, A, user_grad, s, g, beta_0, tilde_m2, ::Val{:TPS}) where {T}
#   D = Descriptor(6, 2)
#   Δ = @vars(D)
#   z = [q[1], p[1], q[2], p[2], q[3], p[3]] .+ Δ
#   h = H_implicit(z, A, s, g, beta_0, tilde_m2)
#   Q = (XI, YI, ZI)
#   P = (PXI, PYI, PZI)
#   c11 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[1]), Q[1])))
#   c12 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[2]), Q[1])))
#   c13 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[3]), Q[1])))
#   c21 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[1]), Q[2])))
#   c22 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[2]), Q[2])))
#   c23 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[3]), Q[2])))
#   c31 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[1]), Q[3])))
#   c32 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[2]), Q[3])))
#   c33 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, P[3]), Q[3])))
#   return SMatrix{3,3,T}(c11, c12, c13, c21, c22, c23, c31, c32, c33)
# end

# @inline function hessian_qp_source(q::SVector{3,T}, p::SVector{3,T}, A, user_grad, s, g, beta_0, tilde_m2, ::Val{:TPS}) where {T}
#   D = Descriptor(6, 2)
#   Δ = @vars(D)
#   z = [q[1], p[1], q[2], p[2], q[3], p[3]] .+ Δ
#   h = H_implicit(z, A, s, g, beta_0, tilde_m2)
#   Q = (XI, YI, ZI)
#   P = (PXI, PYI, PZI)
#   c11 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[1]), P[1])))
#   c12 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[2]), P[1])))
#   c13 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[3]), P[1])))
#   c21 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[1]), P[2])))
#   c22 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[2]), P[2])))
#   c23 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[3]), P[2])))
#   c31 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[1]), P[3])))
#   c32 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[2]), P[3])))
#   c33 = T(GTPSA.scalar(GTPSA.deriv(GTPSA.deriv(h, Q[3]), P[3])))
#   return SMatrix{3,3,T}(c11, c12, c13, c21, c22, c23, c31, c32, c33)
# end
# ==============================================================#

#=

@inline residual_step_q(q, q0, p0, A, s, g, beta_0, tilde_m2, ds2) =
  q - q0 - ds2 * ∇p_H(A, q, p0, s, g, beta_0, tilde_m2)

@inline residual_step_p(p, p1, q1, A, J, s, ds2, g, beta_0, tilde_m2) =
  p - p1 + ds2 * ∇q_H(A, J, q1, p, s, g, beta_0, tilde_m2)


function find_root_q(q0::SVector{3,T}, p0::SVector{3,T}, A, s, ds2, g, beta_0, tilde_m2, ::Val{:DerivativeFree}, abstol, reltol, max_iter) where {T}
  x1, x2, x3 = q0
  r = residual_step_q(SVector{3,T}(x1, x2, x3), q0, p0, A, s, ds2, g, beta_0, tilde_m2)
  r1, r2, r3 = r
  b11 = one(T); b12 = zero(T); b13 = zero(T)
  b21 = zero(T); b22 = one(T); b23 = zero(T)
  b31 = zero(T); b32 = zero(T); b33 = one(T)
  abs_tol = T(abstol)
  rel_tol = T(reltol)
  step_tol = T(1e-32)
  damping = T(0.5)
  for _ in 1:max_iter
    norm_x = sqrt(muladd(x1, x1, muladd(x2, x2, x3 * x3)))
    tol = muladd(norm_x, rel_tol, abs_tol)
    if !(muladd(r1, r1, muladd(r2, r2, r3 * r3)) > tol * tol)
      break
    end
    d1 = muladd(b11, r1, muladd(b12, r2, b13 * r3))
    d2 = muladd(b21, r1, muladd(b22, r2, b23 * r3))
    d3 = muladd(b31, r1, muladd(b32, r2, b33 * r3))
    xt1 = x1 - d1; xt2 = x2 - d2; xt3 = x3 - d3
    rt = residual_step_q(SVector{3,T}(xt1, xt2, xt3), q0, p0, A, s, ds2, g, beta_0, tilde_m2)
    rt1, rt2, rt3 = rt
    res_sq = muladd(r1, r1, muladd(r2, r2, r3 * r3))
    res_trial_sq = muladd(rt1, rt1, muladd(rt2, rt2, rt3 * rt3))
    improved = res_trial_sq < res_sq
    xu1, xu2, xu3 = xt1, xt2, xt3
    ru1, ru2, ru3 = rt1, rt2, rt3
    if !improved
      xr1 = x1 - damping * d1; xr2 = x2 - damping * d2; xr3 = x3 - damping * d3
      rr = residual_step_q(SVector{3,T}(xr1, xr2, xr3), q0, p0, A, s, ds2, g, beta_0, tilde_m2)
      rr1, rr2, rr3 = rr
      xu1 = xr1; xu2 = xr2; xu3 = xr3
      ru1 = rr1; ru2 = rr2; ru3 = rr3
    end
    s1 = xu1 - x1; s2 = xu2 - x2; s3 = xu3 - x3
    y1 = ru1 - r1; y2 = ru2 - r2; y3 = ru3 - r3
    By1 = muladd(b11, y1, muladd(b12, y2, b13 * y3))
    By2 = muladd(b21, y1, muladd(b22, y2, b23 * y3))
    By3 = muladd(b31, y1, muladd(b32, y2, b33 * y3))
    stB1 = muladd(s1, b11, muladd(s2, b21, s3 * b31))
    stB2 = muladd(s1, b12, muladd(s2, b22, s3 * b32))
    stB3 = muladd(s1, b13, muladd(s2, b23, s3 * b33))
    stBy = muladd(stB1, y1, muladd(stB2, y2, stB3 * y3))
    inv_denom = vifelse(abs(stBy) > step_tol, one(T) / stBy, zero(T))
    u1 = s1 - By1; u2 = s2 - By2; u3 = s3 - By3
    b11 = b11 + u1 * stB1 * inv_denom; b12 = b12 + u1 * stB2 * inv_denom; b13 = b13 + u1 * stB3 * inv_denom
    b21 = b21 + u2 * stB1 * inv_denom; b22 = b22 + u2 * stB2 * inv_denom; b23 = b23 + u2 * stB3 * inv_denom
    b31 = b31 + u3 * stB1 * inv_denom; b32 = b32 + u3 * stB2 * inv_denom; b33 = b33 + u3 * stB3 * inv_denom
    x1, x2, x3 = xu1, xu2, xu3
    r1, r2, r3 = ru1, ru2, ru3
  end
  return SVector{3,T}(x1, x2, x3)
end

function find_root_p(p1::SVector{3,T}, q1::SVector{3,T}, A, s, ds2, g, beta_0, tilde_m2, user_grad, ::Val{:DerivativeFree}, abstol, reltol, max_iter) where {T}
  x1, x2, x3 = p1
  r = residual_step_p(SVector{3,T}(x1, x2, x3), p1, q1, A, s, ds2, g, beta_0, tilde_m2, user_grad)
  r1, r2, r3 = r
  b11 = one(T); b12 = zero(T); b13 = zero(T)
  b21 = zero(T); b22 = one(T); b23 = zero(T)
  b31 = zero(T); b32 = zero(T); b33 = one(T)
  abs_tol = T(abstol)
  rel_tol = T(reltol)
  step_tol = T(1e-32)
  damping = T(0.5)
  for _ in 1:max_iter
    norm_x = sqrt(muladd(x1, x1, muladd(x2, x2, x3 * x3)))
    tol = muladd(norm_x, rel_tol, abs_tol)
    if !(muladd(r1, r1, muladd(r2, r2, r3 * r3)) > tol * tol)
      break
    end
    d1 = muladd(b11, r1, muladd(b12, r2, b13 * r3))
    d2 = muladd(b21, r1, muladd(b22, r2, b23 * r3))
    d3 = muladd(b31, r1, muladd(b32, r2, b33 * r3))
    xt1 = x1 - d1; xt2 = x2 - d2; xt3 = x3 - d3
    rt = residual_step_p(SVector{3,T}(xt1, xt2, xt3), p1, q1, A, s, ds2, g, beta_0, tilde_m2, user_grad)
    rt1, rt2, rt3 = rt
    res_sq = muladd(r1, r1, muladd(r2, r2, r3 * r3))
    res_trial_sq = muladd(rt1, rt1, muladd(rt2, rt2, rt3 * rt3))
    improved = res_trial_sq < res_sq
    xu1, xu2, xu3 = xt1, xt2, xt3
    ru1, ru2, ru3 = rt1, rt2, rt3
    if !improved
      xr1 = x1 - damping * d1; xr2 = x2 - damping * d2; xr3 = x3 - damping * d3
      rr = residual_step_p(SVector{3,T}(xr1, xr2, xr3), p1, q1, A, s, ds2, g, beta_0, tilde_m2, user_grad)
      rr1, rr2, rr3 = rr
      xu1 = xr1; xu2 = xr2; xu3 = xr3
      ru1 = rr1; ru2 = rr2; ru3 = rr3
    end
    s1 = xu1 - x1; s2 = xu2 - x2; s3 = xu3 - x3
    y1 = ru1 - r1; y2 = ru2 - r2; y3 = ru3 - r3
    By1 = muladd(b11, y1, muladd(b12, y2, b13 * y3))
    By2 = muladd(b21, y1, muladd(b22, y2, b23 * y3))
    By3 = muladd(b31, y1, muladd(b32, y2, b33 * y3))
    stB1 = muladd(s1, b11, muladd(s2, b21, s3 * b31))
    stB2 = muladd(s1, b12, muladd(s2, b22, s3 * b32))
    stB3 = muladd(s1, b13, muladd(s2, b23, s3 * b33))
    stBy = muladd(stB1, y1, muladd(stB2, y2, stB3 * y3))
    inv_denom = vifelse(abs(stBy) > step_tol, one(T) / stBy, zero(T))
    u1 = s1 - By1; u2 = s2 - By2; u3 = s3 - By3
    b11 = b11 + u1 * stB1 * inv_denom; b12 = b12 + u1 * stB2 * inv_denom; b13 = b13 + u1 * stB3 * inv_denom
    b21 = b21 + u2 * stB1 * inv_denom; b22 = b22 + u2 * stB2 * inv_denom; b23 = b23 + u2 * stB3 * inv_denom
    b31 = b31 + u3 * stB1 * inv_denom; b32 = b32 + u3 * stB2 * inv_denom; b33 = b33 + u3 * stB3 * inv_denom
    x1, x2, x3 = xu1, xu2, xu3
    r1, r2, r3 = ru1, ru2, ru3
  end
  return SVector{3,T}(x1, x2, x3)
end
=#

function find_root_q(A, J, q0::SVector{3,T}, p0::SVector{3,T}, s, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter) where {T}
  x = q0
  abs_tol = T(abstol)
  rel_tol = T(reltol)
  I3 = one(SMatrix{3,3,T})
  for _ in 1:max_iter
    r = x .- q0 .- ds2 * ∇p_H(A, x, p0, s, g, beta_0, tilde_m2)
    tol = muladd(sqrt(norm_sq(x)), rel_tol, abs_tol)
    if !(norm_sq(r) > tol * tol)
      break
    end
    Hp_q = hessian_pq(x, p0, A, J, s, g, beta_0, tilde_m2)
    _J = I3 .- ds2 * Hp_q
    x = x .- solve_3x3(_J, r)
  end
  return SVector(x)
end

function find_root_p(A, J, p1::SVector{3,T}, q1::SVector{3,T}, s, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter) where {T}
  x = p1
  abs_tol = T(abstol)
  rel_tol = T(reltol)
  I3 = one(SMatrix{3,3,T})
  for _ in 1:max_iter
    r = x .- p1 .+ ds2 * ∇q_H(A, J, q1, x, s, g, beta_0, tilde_m2)
    tol = muladd(sqrt(norm_sq(x)), rel_tol, abs_tol)
    if !(norm_sq(r) > tol * tol)
      break
    end
    Hq_p = hessian_qp(q1, x, A, J, s, g, beta_0, tilde_m2)
    _J = I3 .+ ds2 * Hq_p
    x = x .- solve_3x3(_J, r)
  end
  return SVector(x)
end

#=
function H_implicit(z, A, s, g, beta_0, tilde_m2)
  x = z[XI]
  px = z[PXI]
  y = z[YI]
  py = z[PYI]
  ct = z[ZI]
  pz = z[PZI]
  rel_p = 1 + pz
  tilde_E2 = rel_p^2 + tilde_m2
  beta = rel_p / (beta_0 * sqrt(tilde_E2))
  t = ct / (beta * C_LIGHT)
  a = A(x, y, s, t)
  kappa = 1 + g * x
  p_kin_x = px - a[2]
  p_kin_y = py - a[3]
  p_s = sqrt(rel_p^2 - p_kin_x^2 - p_kin_y^2)
  return -kappa * p_s - a[1] - a[4] + sqrt(tilde_E2) / beta_0
end

# THIS REMOVES AN ORDER OF TPSA
# Can be circumvented by providing the vector potential gradient to the solver
function F_implicit(z, A, s, g, beta_0, tilde_m2, ds, forward::Bool)
  h = H_implicit(z, A, s, g, beta_0, tilde_m2)
  sgn = forward ? 1 : -1
  return [z[XI] + sgn * ds * GTPSA.deriv(h, PXI);
    z[PXI] + sgn * ds * GTPSA.deriv(h, XI);
    z[YI] + sgn * ds * GTPSA.deriv(h, PYI);
    z[PYI] + sgn * ds * GTPSA.deriv(h, YI);
    z[ZI] + sgn * ds * GTPSA.deriv(h, PZI);
    z[PZI] + sgn * ds * GTPSA.deriv(h, ZI)]
end
=#

@makekernel fastgtpsa=true function symplectic_step_tpsa!(i, coords::Coords,
  A, J, s, g, beta_0, tilde_m2,
  abstol, reltol, max_iter, ds)
  v = coords.v
  ds2 = ds / 2
  ds4 = ds / 4
  s1 = s[1] + ds4
  s2 = s[1] + 3 * ds4
  T = eltype(v).parameters[1]
  D = eltype(v).parameters[2]
  dz = @vars(D)

  q0_s = SVector{3,T}(GTPSA.scalar(v[i, XI]), GTPSA.scalar(v[i, YI]), GTPSA.scalar(v[i, ZI]))
  p0_s = SVector{3,T}(GTPSA.scalar(v[i, PXI]), GTPSA.scalar(v[i, PYI]), GTPSA.scalar(v[i, PZI]))

  # ======  Part 1 ======

  q1_s = find_root_q(A, J, q0_s, p0_s, s1, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  p1_s = p0_s - ds2 * ∇q_H(A, J, q1_s, p0_s, s1, g, beta_0, tilde_m2)


  p0 = SVector(p0_s[1]+dz[PXI], p0_s[2]+dz[PYI], p0_s[3]+dz[PZI])
  q1 = SVector(q1_s[1]+dz[ XI], q1_s[2]+dz[ YI], q1_s[3]+dz[ ZI])

  ∂qF = q1 - ds2 * ∇p_H(A, q1, p0, s1, g, beta_0, tilde_m2)
  ∂pF = p0 - ds2 * ∇q_H(A, J, q1, p0, s1, g, beta_0, tilde_m2)

  ∂qp_inv = [
    ∂qF[1],∂pF[1],
    ∂qF[2],∂pF[2],
    ∂qF[3],∂pF[3]
  ]

  ∂qp = zero(∂qp_inv)
  
  GTPSA.pminv!(6, ∂qp_inv, 6, ∂qp, [1, 0, 1, 0, 1, 0])

  v[i, :] = [q1_s[1], p1_s[1], q1_s[2], p1_s[2], q1_s[3], p1_s[3]] .+ (∂qp ∘ cutord.(v[i, :], 0))

  # ======  Part 2 ======

  p2_s = find_root_p(A, J, p1_s, q1_s, s2, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  q2_s = q1_s + ds2 * ∇p_H(A, q1_s, p2_s, s2, g, beta_0, tilde_m2)

  p2 = SVector(p2_s[1]+dz[PXI], p2_s[2]+dz[PYI], p2_s[3]+dz[PZI])
  q1 = SVector(q1_s[1]+dz[ XI], q1_s[2]+dz[ YI], q1_s[3]+dz[ ZI])

  ∂qF = q1 + ds2 * ∇p_H(A, q1, p2, s2, g, beta_0, tilde_m2)
  ∂pF = p2 + ds2 * ∇q_H(A, J, q1, p2, s2, g, beta_0, tilde_m2)

  ∂qp_inv = [
    ∂qF[1],∂pF[1],
    ∂qF[2],∂pF[2],
    ∂qF[3],∂pF[3]
  ]

  #∂qp = zero(∂qp_inv)

  GTPSA.pminv!(6, ∂qp_inv, 6, ∂qp, [0, 1, 0, 1, 0, 1])

  v[i, :] = [q2_s[1], p2_s[1], q2_s[2], p2_s[2], q2_s[3], p2_s[3]] .+ (∂qp ∘ cutord.(v[i, :], 0))


  s[1] += ds
end

@makekernel function symplectic_step!(i, coords::Coords,
  A, J, s, a, g, beta_0, tilde_m2,
  abstol, reltol, max_iter, ds)
  v = coords.v
  T = eltype(v)
  q0 = SVector{3,T}(v[i, XI], v[i, YI], v[i, ZI])
  p0 = SVector{3,T}(v[i, PXI], v[i, PYI], v[i, PZI])
  ds2 = ds / 2
  ds4 = ds / 4
  s1 = s[1] + ds4
  s2 = s[1] + 3 * ds4

  #s[1] += ds4

  q1 = find_root_q(A, J, q0, p0, s1, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  p1 = p0 .- ds2 * ∇q_H(A, J, q1, p0, s1, g, beta_0, tilde_m2)

  #s[1] += ds4

  #s[1] += ds4

  p2 = find_root_p(A, J, p1, q1, s2, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  q2 = q1 .+ ds2 * ∇p_H(A, q1, p2, s2, g, beta_0, tilde_m2)

  #s[1] += ds4

  s[1] += ds

  alive = (coords.state[i] == STATE_ALIVE)
  v[i, XI] = vifelse(alive, q2[1], v[i, XI])
  v[i, PXI] = vifelse(alive, p2[1], v[i, PXI])
  v[i, YI] = vifelse(alive, q2[2], v[i, YI])
  v[i, PYI] = vifelse(alive, p2[2], v[i, PYI])
  v[i, ZI] = vifelse(alive, q2[3], v[i, ZI])
  v[i, PZI] = vifelse(alive, p2[3], v[i, PZI])

  if !isnothing(coords.q)
    rel_p0 = 1 + p0[3]; rel_p2 = 1 + p2[3]
    tilde_E0_sq = muladd(rel_p0, rel_p0, tilde_m2)
    tilde_E2_sq = muladd(rel_p2, rel_p2, tilde_m2)
  
    beta0_inv = beta_0 * sqrt(tilde_E0_sq)
    beta2_inv = beta_0 * sqrt(tilde_E2_sq)
    beta0 = rel_p0 / beta0_inv
    beta2 = rel_p2 / beta2_inv
    t0 = q0[3] / (beta0 * C_LIGHT)
    t2 = q2[3] / (beta2 * C_LIGHT)

    A0 = A(q0[1],q0[2],s1-ds4,t0);    A2 = A(q2[1],q2[2],s2+ds4,t2);
    J0 = J(q0[1],q0[2],s1-ds4,t0);    J2 = J(q2[1],q2[2],s2+ds4,t2);
    
    E0 = SVector(J0[1,1], J0[1,2], J0[1,3])
    B0 = SVector(J0[4,2]-J0[3,3], J0[2,3]-J0[4,1], J0[3,1]-J0[2,2])

    E2 = SVector(J2[1,1], J2[1,2], J2[1,3])
    B2 = SVector(J2[4,2]-J2[3,3], J2[2,3]-J2[4,1], J2[3,1]-J2[2,2])
    
    Ω0 = omega_field(i, coords, a, g, sqrt(tilde_m2), A0[1], A0[2], E0, B0, ds / 2)
    Ω2 = omega_field(i, coords, a, g, sqrt(tilde_m2), A2[1], A2[2], E2, B2, ds / 2)

    coords.q[i,Q0], coords.q[i,QX], coords.q[i,QY], coords.q[i,QZ] = quat_mul(
      expq( Ω0 .+ Ω2, alive ), 
      coords.q[i,Q0], coords.q[i,QX], coords.q[i,QY], coords.q[i,QZ]
      )
  end
end
