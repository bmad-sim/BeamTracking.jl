@makekernel fastgtpsa=true function track_alignment_straight!(i, coords::Coords, beta_0, gamsqr_0, tilde_m,
       entering, x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)

  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  L2 = 0.5 * L * ele_orient

  if entering
    v[i,XI] = vifelse(alive, v[i,XI] - x_off, v[i,XI])
    v[i,YI] = vifelse(alive, v[i,YI] - y_off, v[i,YI])
 
    q = inv_rot_quaternion(x_rot, y_rot, tilt)
    dz = -L2 - z_off    # Z distance from element center
    dz_new = coord_rotation!(i, coords, q, dz)

    ExactTracking.exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, dz_new + L2)

  else
    q = rot_quaternion(x_rot, y_rot, tilt)
    dz = L2 + z_off
    dz_new = coord_rotation!(i, coords, q, dz)

    v[i,XI] = vifelse(alive, v[i,XI] + x_off, v[i,XI])
    v[i,YI] = vifelse(alive, v[i,YI] + y_off, v[i,YI])

    ExactTracking.exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, dz_new - L2)
  end
end
 
#

@makekernel fastgtpsa=true function track_alignment_bend!(i, coords::Coords, beta_0, gamsqr_0, tilde_m,
                    entering, x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)
  v = coords.v

end
 
