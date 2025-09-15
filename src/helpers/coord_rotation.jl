"""
    coord_rotation!(i, coords::Coords, q_inv, z_0) 
Rotates both `(x, y)` and `(px, py, pz)` phase space coordinates.

`z_0` is the longitudinal offset from the center of rotation. 

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

    w31 =     2*(q_inv[QX]*q_inv[QZ] + q_inv[QY]*q_inv[Q0])
    w32 =     2*(q_inv[QY]*q_inv[QZ] + q_inv[QX]*q_inv[Q0])
    w33 = 1 - 2*(q_inv[QX]*q_inv[QX] + q_inv[QY]*q_inv[QY])

    x_0 = v[i,XI]
    y_0 = v[i,YI]
    new_x = w11*x_0 + w12*y_0 - w13*z_0
    new_y = w21*x_0 + w22*y_0 - w23*z_0
    new_z = w31*x_0 + w32*y_0 - w33*z_0

    v[i,XI] = vifelse(alive, new_x, x_0)
    v[i,YI] = vifelse(alive, new_y, y_0)
    z_out = vifelse(alive, new_z, z_0)

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
