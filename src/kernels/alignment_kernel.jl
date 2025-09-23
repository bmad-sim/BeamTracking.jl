@makekernel fastgtpsa=true function track_alignment_straight_entering!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  L2 = 0.5 * L * ele_orient

  v[i,XI] = vifelse(alive, v[i,XI] - x_off, v[i,XI])
  v[i,YI] = vifelse(alive, v[i,YI] - y_off, v[i,YI])
 
  q = inv_rot_quaternion(x_rot, y_rot, tilt)
  dz_new = coord_rotation!(i, coords, q, -L2 - z_off)
  alignment_drift!(i, coords, -dz_new - L2)
end
 
#

@makekernel fastgtpsa=true function track_alignment_straight_exiting!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  L2 = 0.5 * L * ele_orient

  q = rot_quaternion(x_rot, y_rot, tilt)
  dz_new = coord_rotation!(i, coords, q, L2)

  v[i,XI] = vifelse(alive, v[i,XI] + x_off, v[i,XI])
  v[i,YI] = vifelse(alive, v[i,YI] + y_off, v[i,YI])

  alignment_drift!(i, coords, L2 - z_off - dz_new)
end
 
#

@makekernel fastgtpsa=true function track_alignment_bend_entering!(i, coords::Coords,
                          x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  # Transform to coord system with origin at element chord midpoint.

  z1 = coord_bend_transform!(i, coords, 0.5*L, g_ref, tilt_ref, 0.0)
  dr = -L^2 * g_ref * one_cos_norm(tilt) * [cos(tilt_ref), sin(tilt_ref), 0.0]
  z1 = coord_translation!(i, coords, dr, z1)

  # Misalignment transform

  v[i,XI] = vifelse(alive, v[i,XI] - x_off, v[i,XI])
  v[i,YI] = vifelse(alive, v[i,YI] - y_off, v[i,YI])
 
  q = inv_rot_quaternion(x_rot, y_rot, tilt)
  z1 = coord_rotation!(i, coords, q, z1)

  # Transform to body coord system with origin at entrance.

  q = inv_rot_quaternion(0.0, 0.0, tilt_ref)
  z1 = coord_rotation!(i, coords, q, z1)

  dr = L^2 * g_ref * one_cos_norm(tilt) * [1.0, 0.0, 0.0]
  z1 = coord_translation!(i, coords, dr, z1)
  z1 = coord_bend_transform!(i, coords, -0.5*L, g_ref, 0.0, z1)

  # And drift to entrance

  alignment_drift!(i, coords, -z1)
end
 
#

@makekernel fastgtpsa=true function track_alignment_bend_exiting!(i, coords::Coords,
                          x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  L2 = 0.5 * L * ele_orient


end
 