"""
This function computes sin(sqrt(x))/sqrt(x) and cos(sqrt(x)), which are both 
necessary for exponentiating a rotation vector into a quaternion.
"""
function sincos_quaternion(x)
  threshold = 7.3e-8 # sqrt(24*eps(Float64))
  sq = sqrt(x)
  s, c = sincos(sq)
  s = s/sq
  s_out = vifelse(x > threshold, s, 1-x/6)
  c_out = vifelse(x > threshold, c, 1-x/2)
  return s_out, c_out
end


"""
This function computes sin(sqrt(x))/sqrt(x) and cos(sqrt(x)), which are both 
necessary for exponentiating a rotation vector into a quaternion.
"""
function sincos_quaternion(x::TPS{T}) where {T}
  ε = eps(T)
  N_max = 100
  N = 1
  conv_sin = false
  conv_cos = false
  y = one(x)
  prev_sin = one(x)
  prev_cos = one(x)
  result_sin = one(x)
  result_cos = one(x)
  #sq = one(x)
  # Using FastGTPSA! for the following makes other kernels run out of temps
  @FastGTPSA begin
    if x < ε #0.1
      while !(conv_sin && conv_cos) && N < N_max
        y = -y*x/((2*N)*(2*N - 1))
        result_sin = prev_sin + y/(2*N + 1)
        result_cos = prev_cos + y
        N += 1
        if normTPS(result_sin - prev_sin) < ε
          conv_sin = true
        end
        if normTPS(result_cos - prev_cos) < ε
          conv_cos = true
        end
        prev_sin = result_sin
        prev_cos = result_cos
      end
    else
      sq = sqrt(x)
      result_sin, result_cos = sincos(sq)
      result_sin = result_sin/sq
    end
  end
  if N == N_max
    @warn "sincos_quaternion convergence not reached in $N_max iterations"
  end
  return result_sin, result_cos
end


"""
This function computes exp(-i/2 v⋅σ) as a quaternion, where σ is the 
vector of Pauli matrices. If compute is false, it returns the identity quaternion.
"""
function expq(v, compute)
  n2 = @FastGTPSA (v[1]*v[1] + v[2]*v[2] + v[3]*v[3])/4
  n2_0 = zero(n2)
  s, c = sincos_quaternion(vifelse(compute, n2, n2_0))
  s = vifelse(compute, s, n2_0)
  return (c, s*v[1]/2, s*v[2]/2, s*v[3]/2)
end


"""
    quat_mul(q1, q20, q2x, q2y, q2z) -> q3 = q1*q2

Quaternion multiplication`q1 * q2` where `q2 = [q20, q2x, q2y, q2z]`.
This form of `quat_mul` is used when the quaternions are particle (spin) coordinates and is needed
with SIMD-parallelized tracking.
"""
function quat_mul(q1, q20, q2x, q2y, q2z)
  a1, b1, c1, d1 = q1[Q0], q1[QX], q1[QY], q1[QZ]
  a2, b2, c2, d2 = q20, q2x, q2y, q2z
  @FastGTPSA begin
    a3 = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b3 = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c3 = a1*c2 - b1*d2 + c1*a2 + d1*b2
    d3 = a1*d2 + b1*c2 - c1*b2 + d1*a2
  end
  return (a3, b3, c3, d3)
end


"""
    quat_mul(q1, q2) -> q3 = q1*q2

Quaternion multiplication`q1 * q2`.

Also see `quat_mul(q1, q20, q2x, q2y, q2z)` which iss the form of `quat_mul` that is used when 
the quaternions are particle (spin) coordinates and is needed with SIMD-parallelized tracking.
"""
@inline quat_mul(q1, q2) = quat_mul(q1, q2[Q0], q2[QX], q2[QY], q2[QZ])


@inline quat_inv(q) = (q[Q0], -q[QX], -q[QY], -q[QZ])



"""
    function quat_rotate(r, q) -> rotated_r

Rotates vector `r` using quaternion `q`.
"""
@inline function quat_rotate(r, q)

  w11 = 1 - 2*(q[QY]*q[QY] + q[QZ]*q[QZ])
  w12 =     2*(q[QX]*q[QY] - q[QZ]*q[Q0])
  w13 =     2*(q[QX]*q[QZ] + q[QY]*q[Q0])

  w21 =     2*(q[QX]*q[QY] + q[QZ]*q[Q0])
  w22 = 1 - 2*(q[QX]*q[QX] + q[QZ]*q[QZ])
  w23 =     2*(q[QY]*q[QZ] - q[QX]*q[Q0])

  w31 =     2*(q[QX]*q[QZ] - q[QY]*q[Q0])
  w32 =     2*(q[QY]*q[QZ] + q[QX]*q[Q0])
  w33 = 1 - 2*(q[QX]*q[QX] + q[QY]*q[QY])

  return (w11*r[1] + w12*r[2] + w13*r[3],
          w21*r[1] + w22*r[2] + w23*r[3],
          w31*r[1] + w32*r[2] + w33*r[3])
end


"""
    rot_quat(axis, angle) -> q

Calculates rotation quaternion from axis and angle.
It is assumed that the axis is properly normalized.

## Arguments
- `axis`      Three vector axis.
- `angle`     Rotation angle.

# Returns
- `q`    quaternion 4-vector.
"""
function rot_quat(axis, angle)
  s = sin(0.5*angle)
  return (cos(0.5*angle), s*axis[1], s*axis[2], s*axis[3])
end


"""
  rot_quaternion(x_rot, y_rot, z_rot)

Constructs a rotation quaternion based on the given Bryan-Tait angles.

Bmad/SciBmad follows the MAD convention of applying z, x, y rotations in that order.

The inverse quaternion reverses the order of operations and their signs.


Arguments:
- `x_rot::Number`: Rotation angle around the x-axis.
- `y_rot::Number`: Rotation angle around the y-axis.
- `z_rot::Number`: Rotation angle around the z-axis.

"""
function rot_quaternion(x_rot, y_rot, z_rot)
  qz = SA[cos(z_rot/2), 0, 0, sin(z_rot/2)]
  qx = SA[cos(x_rot/2), sin(x_rot/2), 0, 0]
  qy = SA[cos(y_rot/2), 0, sin(y_rot/2), 0]
  q = quat_mul(qx, qz[Q0], qz[QX], qz[QY], qz[QZ])
  q = quat_mul(qy, q[Q0], q[QX], q[QY], q[QZ])
  return SA[q[Q0], q[QX], q[QY], q[QZ]]
end

# Inverse rotation quaternion
function inv_rot_quaternion(x_rot, y_rot, z_rot)
  qz = SA[cos(z_rot/2), 0, 0, -sin(z_rot/2)]
  qx = SA[cos(x_rot/2), -sin(x_rot/2), 0, 0]
  qy = SA[cos(y_rot/2), 0, -sin(y_rot/2), 0]
  q = quat_mul(qx, qy[Q0], qy[QX], qy[QY], qy[QZ])
  q = quat_mul(qz, q[Q0], q[QX], q[QY], q[QZ])
  return SA[q[Q0], q[QX], q[QY], q[QZ]]
end
