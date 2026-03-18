"""
    norm_sq(v)

Return the Euclidean norm squared of a three-vector `v`.
"""
@inline norm_sq(v::SVector{3,T}) where {T} = muladd(v[1], v[1], muladd(v[2], v[2], v[3] * v[3]))

"""
    ∇q_H(A, Jac, q, p, s, g, beta_0, tilde_m2)

Return the gradient of the implicit Hamiltonian with respect to
the canonical coordinates `q = (x, y, z)` at the phase-space point
`(q, p)`. The four-potential `A` and its Jacobian `Jac` are
evaluated at the longitudinal position `s` and the corresponding
time of flight.

## Arguments
- `A`:         vector potential callback returning `(ϕ, Ax, Ay, As)`
- `Jac`:       Jacobian callback for `A`
- `q`:         canonical coordinates `(x, y, z)`
- `p`:         canonical momenta `(p_x, p_y, p_z)`
- `s`:         longitudinal integration position
- `g`:         reference curvature
- `beta_0`:    relativistically-normalized speed
- `tilde_m2`:  squared, normalized rest mass `tilde_m2 == (mc/p)^2`
"""
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

  dH_dqx = -_J[4, 1] - scale * muladd(p_kin_x, _J[2, 1], p_kin_y * _J[3, 1]) - g * p_s
  dH_dqy = -_J[4, 2] - scale * muladd(p_kin_x, _J[2, 2], p_kin_y * _J[3, 2])
  dH_dqz = -_J[4, 3] - scale * muladd(p_kin_x, _J[2, 3], p_kin_y * _J[3, 3])

  return SVector{3,T}(dH_dqx, dH_dqy, dH_dqz)
end

"""
    ∇p_H(A, q, p, s, g, beta_0, tilde_m2)

Return the gradient of the implicit Hamiltonian with respect to
the canonical momenta `p = (p_x, p_y, p_z)` at the phase-space point
`(q, p)`.

## Arguments
- `A`:         vector potential callback returning `(ϕ, Ax, Ay, As)`
- `q`:         canonical coordinates `(x, y, z)`
- `p`:         canonical momenta `(p_x, p_y, p_z)`
- `s`:         longitudinal integration position
- `g`:         reference curvature
- `beta_0`:    relativistically-normalized speed
- `tilde_m2`:  squared, normalized rest mass `tilde_m2 == (mc/p)^2`
"""
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

"""
    solve_3x3(A, b) == A / b

Solve the dense linear system `A * x = b` for a 3-vector `x`.
"""
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


"""
    mixed_hessian_pq(q, p, A, Jac, s, g, beta_0, tilde_m2)

Return the mixed Hessian of the implicit Hamiltonian evaluated at
`(q, p)`, with momentum components indexing rows and coordinate
components indexing columns:

`Hpq[i, j] = ∂²H / (∂pᵢ ∂qⱼ)`.

Equivalently, row 1/2/3 corresponds to `(p_x, p_y, p_z)` and
column 1/2/3 corresponds to `(x, y, z)`.
"""
@inline function mixed_hessian_pq(q::SVector{3,T}, p, A, Jac, s, g, beta_0, tilde_m2) where {T}
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


"""
    find_root_q(A, J, q0, p0, s, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)

Solve the implicit position half-step

`q1 = q0 + ds2 * ∂H/∂p(q1, p0)`

with Newton iteration.
"""
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
    Hp_q = mixed_hessian_pq(x, p0, A, J, s, g, beta_0, tilde_m2)
    _J = I3 .- ds2 * Hp_q
    x = x .- solve_3x3(_J, r)
  end
  return SVector(x)
end

"""
    find_root_p(A, J, p1, q1, s, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)

Solve the implicit momentum half-step

`p2 = p1 - ds2 * ∂H/∂q(q1, p2)`

with Newton iteration.
"""
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
    # `mixed_hessian_pq` stores p-components on rows and q-components on
    # columns, so transpose it here to get q-rows / p-columns.
    Hq_p = SMatrix{3,3,T,9}(transpose(mixed_hessian_pq(q1, x, A, J, s, g, beta_0, tilde_m2)))
    _J = I3 .+ ds2 * Hq_p
    x = x .- solve_3x3(_J, r)
  end
  return SVector(x)
end


"""
symplectic_step!()

Track a particle through one implicit symplectic step using a
position-implicit half-step followed by a momentum-implicit
half-step. When spin coordinates are present, the routine also
applies a second-order spin rotation using the endpoint fields.

The orbital step is written as the composition of two canonical maps.
The first half-step is generated by the type-3 generating function
`F₃(p₀, q₁) = -q₁⋅p₀ + (ds/2) H(q₁, p₀)`, which yields the implicit
solve for `q₁` and the matching update for `p₁`. The second half-step
is generated by the type-2 generating function
`F₂(q₁, p₂) = q₁⋅p₂ + (ds/2) H(q₁, p₂)`.

These two substeps are adjoints of each other, so the symmetric
composition `F₃(ds/2) ∘ F₂(ds/2)`, evaluated at `s + ds/4` and
`s + 3ds/4`, is self-adjoint. That symmetry is what makes the full
orbital map time-reversible and second-order accurate.

Arguments
—————————
A:          vector potential callback returning `(ϕ, Ax, Ay, As)`
J:          Jacobian callback for `A`
s:          singleton container for the initial longitudinal position
a:          gyromagnetic anomaly
g:          reference curvature
beta_0:     relativistically-normalized speed
tilde_m2:   squared, normalized rest mass `tilde_m2 == (mc/p)^2`
abstol:     absolute Newton tolerance
reltol:     relative Newton tolerance
max_iter:   maximum Newton iterations per half-step
ds:         integrator step length
"""
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

  # First half-step: type-3 generating function
  #   F₃(p₀, q₁) = -q₁⋅p₀ + (ds/2) H(q₁, p₀).
  q1 = find_root_q(A, J, q0, p0, s1, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  p1 = p0 .- ds2 * ∇q_H(A, J, q1, p0, s1, g, beta_0, tilde_m2)


  # Second half-step: type-2 generating function
  #   F₂(q₁, p₂) = q₁⋅p₂ + (ds/2) H(q₁, p₂).
  p2 = find_root_p(A, J, p1, q1, s2, ds2, g, beta_0, tilde_m2, abstol, reltol, max_iter)
  q2 = q1 .+ ds2 * ∇p_H(A, q1, p2, s2, g, beta_0, tilde_m2)

  s[1] += ds

  rel_p2 = 1 + p2[3]
  tilde_E2_sq = muladd(rel_p2, rel_p2, tilde_m2)
  beta2_inv = beta_0 * sqrt(tilde_E2_sq)
  beta2 = rel_p2 / beta2_inv
  t2 = q2[3] / (beta2 * C_LIGHT)
  A2 = A(q2[1], q2[2], s2 + ds4, t2)


  if !isnothing(coords.q)
    rel_p0 = 1 + p0[3]
    tilde_E0_sq = muladd(rel_p0, rel_p0, tilde_m2)
    beta0_inv = beta_0 * sqrt(tilde_E0_sq)
    beta0 = rel_p0 / beta0_inv
    t0 = q0[3] / (beta0 * C_LIGHT)
    tilde_m = sqrt(tilde_m2)

    A0 = A(q0[1], q0[2], s1 - ds4, t0)
    J0 = J(q0[1], q0[2], s1 - ds4, t0)
    J2 = J(q2[1], q2[2], s2 + ds4, t2)
    
    E0 = SVector(J0[1,1], J0[1,2], J0[1,3])
    B0 = SVector(J0[4,2]-J0[3,3], J0[2,3]-J0[4,1], J0[3,1]-J0[2,2])

    E2 = SVector(J2[1,1], J2[1,2], J2[1,3])
    B2 = SVector(J2[4,2]-J2[3,3], J2[2,3]-J2[4,1], J2[3,1]-J2[2,2])
    
    Ω0, _ = omega_field_core(q0[1], p0[1], p0[2], p0[3], a, g, tilde_m, A0[2], A0[3], E0, B0, ds / 2)
    Ω2, good_momenta2 = omega_field_core(q2[1], p2[1], p2[2], p2[3], a, g, tilde_m, A2[2], A2[3], E2, B2, ds / 2)
    Ωsum = SVector(Ω0[1] + Ω2[1], Ω0[2] + Ω2[2], Ω0[3] + Ω2[3])
  else
    pkin_x2 = p2[1] - A2[2]
    pkin_y2 = p2[2] - A2[3]
    pl2_2 = muladd(rel_p2, rel_p2, -(muladd(pkin_x2, pkin_x2, pkin_y2 * pkin_y2)))
    good_momenta2 = (pl2_2 > zero(pl2_2))
  end

  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta2 & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)

  v[i, XI] = vifelse(alive, q2[1], v[i, XI])
  v[i, PXI] = vifelse(alive, p2[1], v[i, PXI])
  v[i, YI] = vifelse(alive, q2[2], v[i, YI])
  v[i, PYI] = vifelse(alive, p2[2], v[i, PYI])
  v[i, ZI] = vifelse(alive, q2[3], v[i, ZI])
  v[i, PZI] = vifelse(alive, p2[3], v[i, PZI])

  if !isnothing(coords.q)
    coords.q[i,Q0], coords.q[i,QX], coords.q[i,QY], coords.q[i,QZ] = quat_mul(
      expq(Ωsum, alive), 
      coords.q[i,Q0], coords.q[i,QX], coords.q[i,QY], coords.q[i,QZ]
      )
  end
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
# This is circumvented by providing the vector potential gradient to the solver
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

"""
symplectic_step_tpsa!()

Track a TPSA coordinate map through one implicit symplectic step
using a position-implicit half-step followed by a momentum-implicit
half-step. The scalar fixed point is solved first, after which the
polynomial map is reconstructed by implicit inversion.

As in the scalar tracker, the first half-step is generated by the
type-3 generating function
`F₃(p₀, q₁) = -q₁⋅p₀ + (ds/2) H(q₁, p₀)`, while the second half-step
is generated by the type-2 generating function
`F₂(q₁, p₂) = q₁⋅p₂ + (ds/2) H(q₁, p₂)`.

Because these two substeps are adjoints, their symmetric composition
defines a self-adjoint TPSA map. That inherited symmetry makes the
resulting map time-reversible and second-order accurate.

Arguments
—————————
A:          vector potential callback returning `(ϕ, Ax, Ay, As)`
J:          Jacobian callback for `A`
s:          singleton container for the initial longitudinal position
a:          gyromagnetic anomaly
g:          reference curvature
beta_0:     relativistically-normalized speed
tilde_m2:   squared, normalized rest mass `tilde_m2 == (mc/p)^2`
abstol:     absolute Newton tolerance
reltol:     relative Newton tolerance
max_iter:   maximum Newton iterations per half-step
ds:         integrator step length
"""
@makekernel fastgtpsa=true function symplectic_step_tpsa!(i, coords::Coords,
  A, J, s, a, g, beta_0, tilde_m2,
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

  # First half-step: type-3 generating function.
  #   F₃(p₀, q₁) = -q₁⋅p₀ + (ds/2) H(q₁, p₀).
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

  # Second half-step: type-2 generating function.
  #   F₂(q₁, p₂) = q₁⋅p₂ + (ds/2) H(q₁, p₂).
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
