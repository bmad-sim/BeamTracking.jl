#---------------------------------------------------------------------------------------------------
# coord_transform_concatenation(dr1, q1, dr2, q2) -> dr, q

"""
    coord_transform_concatenation(dr1, q1, dr2, q2) -> dr, q

Returns the composite transform for the coordinate transform (`dr1`, `q1`) followed by (`dr2`, `q2`).

## Output
- `dr`    Composite coordinate origin shift.
- `q`     Composite Quaternion rotation.
"""
function coord_transform_concatenation(dr1, q1, dr2, q2)
  return quat_rot!(dr1, q2) + dr2, quat_mul(q1, q2)
end

#---------------------------------------------------------------------------------------------------
# coord_bend_transform

"""
    function coord_bend_transform(ds, g_ref, tilt_ref) -> dr, q

Return the transform to the particle phase-space coordinates when the coordinate system is 
transformed along a bend arc a distance `ds`.

## Arguments
- `ds`:       Distance along the bend arc the coordinate system is transformed.
- `g_ref`     Reference `g = 1/radius`.
- `tilt_ref`  Reference tilt.

## Output
- `dr`    Coordinate origin shift.
- `q`     Quaternion rotation.
"""
@inline function coord_bend_transform(ds, g_ref, tilt_ref)
  @FastGTPSA begin @inbounds begin 
    ang = ds * g_ref
    sinc_ang = sincu(0.5*ang)
    axis = [sin(tilt_ref), cos(tilt_ref), 0.0]

    q = rot_q(axis, ds * g_ref)
    dr = [0.5 * ds * ang * sinc_ang, 0.0, -ds * sincu(ang)]

    return dr, q
  end end
end

#---------------------------------------------------------------------------------------------------
"""
    function coord_translation!(i, coords::Coords, dr, z0) -> z1

Transform phase-space particle coordinates when the coordinate system is translated by the vector `dr`.

## Arguments
- `ds`:       Distance along the bend arc the coordinate system is transformed.
- `g_ref`     Reference `g = 1/radius`.
- `tilt_ref`  Reference tilt.
- `z0`       The particle longitudinal distance from the center of rotation. 

## Output
- `z1`        The particle longitudinal distance from the center of rotation after rotation. 

"""
@inline function coord_translation!(i, coords::Coords, dr, z0)
  @FastGTPSA begin @inbounds begin 
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)

    x = v[i,XI] - dr[1]
    y = v[i,YI] - dr[2]
    z = z0    - dr[3]

    v[i,XI] = vifelse(alive, x, v[i,XI])
    v[i,YI] = vifelse(alive, y, v[i,YI])
    return z
  end end
end

#---------------------------------------------------------------------------------------------------

"""
    coord_rotation!(i, coords::Coords, q_inv, z_0) -> z_1

Particle coordinate rotation. This rotates both `(x, y)` and `(px, py, pz)` phase space coordinates.

## Arguments
- `q_inv`   The particle quaternion rotation. This is the inverse of the rotation of the coordinate system.
- `z_0`     The longitudinal particle distance from the center of rotation. 

## Output
- `z1`        The particle longitudinal distance from the center of rotation after rotation. 

Returned is the new longitudinal offset from the center of rotation.
"""

@inline function coord_rotation!(i, coords::Coords, q_inv, z_0) 
  @FastGTPSA begin @inbounds begin
    v = coords.v
    rel_p = 1 + v[i,PZI]
    ps_02 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
    good_momenta = (ps_02 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps_02_1 = one(ps_02)
    ps_0 = sqrt(vifelse(good_momenta, ps_02, ps_02_1))

    w11 = 1 - 2*(q_inv[QY]*q_inv[QY] + q_inv[QZ]*q_inv[QZ])
    w12 =     2*(q_inv[QX]*q_inv[QY] - q_inv[QZ]*q_inv[Q0])
    w13 =     2*(q_inv[QX]*q_inv[QZ] + q_inv[QY]*q_inv[Q0])

    w21 =     2*(q_inv[QX]*q_inv[QY] + q_inv[QZ]*q_inv[Q0])
    w22 = 1 - 2*(q_inv[QX]*q_inv[QX] + q_inv[QZ]*q_inv[QZ])
    w23 =     2*(q_inv[QY]*q_inv[QZ] - q_inv[QX]*q_inv[Q0])

    w31 =     2*(q_inv[QX]*q_inv[QZ] - q_inv[QY]*q_inv[Q0])
    w32 =     2*(q_inv[QY]*q_inv[QZ] + q_inv[QX]*q_inv[Q0])
    w33 = 1 - 2*(q_inv[QX]*q_inv[QX] + q_inv[QY]*q_inv[QY])

    x_0 = v[i,XI]
    y_0 = v[i,YI]
    new_x = w11*x_0 + w12*y_0 + w13*z_0
    new_y = w21*x_0 + w22*y_0 + w23*z_0
    new_z = w31*x_0 + w32*y_0 + w33*z_0

    v[i,XI] = vifelse(alive, new_x, x_0)
    v[i,YI] = vifelse(alive, new_y, y_0)
    zero_z = zero(new_z)
    z_out = vifelse(alive, new_z, zero_z)

    px_0 = v[i,PXI]
    py_0 = v[i,PYI]
    new_px = w11*px_0 + w12*py_0 + w13*ps_0
    new_py = w21*px_0 + w22*py_0 + w23*ps_0
    v[i,PXI] = vifelse(alive, new_px, px_0)
    v[i,PYI] = vifelse(alive, new_py, py_0)
  
    q1 = coords.q
    if !isnothing(q1)
      q = quat_mul(q_inv, q1[i,Q0], q1[i,QX], q1[i,QY], q1[i,QZ])
      q0 = vifelse(alive, q[Q0], q1[i,Q0])
      qx = vifelse(alive, q[QX], q1[i,QX])
      qy = vifelse(alive, q[QY], q1[i,QY])
      qz = vifelse(alive, q[QZ], q1[i,QZ])
      q1[i,Q0], q1[i,QX], q1[i,QY], q1[i,QZ] = q0, qx, qy, qz
    end
  
    return z_out
  end end
end
