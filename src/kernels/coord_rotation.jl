@makekernel fastgtpsa=true function coord_rotation!(i, coords::Coords, winv, dz) 
  v = coords.v
  rel_p = 1 + v[i,PZI]
  ps_02 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  good_momenta = (ps_02 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_02_1 = one(ps_02)
  ps_0 = sqrt(vifelse(good_momenta, ps_02, ps_02_1))

  w11 = 1 - 2*(winv[QY]*winv[QY] + winv[QZ]*winv[QZ])
  w12 = 2*(winv[QX]*winv[QY] - winv[QZ]*winv[Q0])
  w13 = 2*(winv[QX]*winv[QZ] + winv[QY]*winv[Q0])
  w21 = 2*(winv[QX]*winv[QY] + winv[QZ]*winv[Q0])
  w22 = 1 - 2*(winv[QX]*winv[QX]+winv[QZ]*winv[QZ])
  w23 = 2*(winv[QY]*winv[QZ] - winv[QX]*winv[Q0])

  x_0 = v[i,XI]
  y_0 = v[i,YI]
  new_x = w11*x_0 + w12*y_0 - w13*dz
  new_y = w21*x_0 + w22*y_0 - w23*dz
  v[i,XI] = vifelse(alive, new_x, x_0)
  v[i,YI] = vifelse(alive, new_y, y_0)

  px_0 = v[i,PXI]
  py_0 = v[i,PYI]
  new_px = w11*px_0 + w12*py_0 + w13*ps_0
  new_py = w21*px_0 + w22*py_0 + w23*ps_0
  v[i,PXI] = vifelse(alive, new_px, px_0)
  v[i,PYI] = vifelse(alive, new_py, py_0)

  q1 = coords.q
  if !isnothing(q1)
    q = quat_mul(winv, q1[i,Q0], q1[i,QX], q1[i,QY], q1[i,QZ])
    q0 = vifelse(alive, q[Q0], q1[i,Q0])
    qx = vifelse(alive, q[QX], q1[i,QX])
    qy = vifelse(alive, q[QY], q1[i,QY])
    qz = vifelse(alive, q[QZ], q1[i,QZ])
    q1[i,Q0], q1[i,QX], q1[i,QY], q1[i,QZ] = q0, qx, qy, qz
  end
end
