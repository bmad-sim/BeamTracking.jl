# entering = 1 -> entering element.
# entering = -1 -> exiting element.

@makekernel fastgtpsa=true function track_alignment_straight!(i, coords::Coords, entering,
          x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, tilde_m, beta_0, L)
  v = coords.v
  L2 = 0.5 * L * ele_orient * entering
  q = rot_quaternion(x_rot, y_rot, tilt)
  r = (v[i,XI] - x_off, v[i,YI] - y_off, -z_off-L2)

  pz2 = (1+v[i,PZI])^2 - v[i,PXI]^2, - v[i,PYI]
  coords.state[i] = vifelse(pz2 < 0, State.Lost, coords.state[i])
  alive = vifelse(coords.state[i]==STATE_ALIVE, 1, 0)
  p = [v[i,PXI], v[i,PYI], sqrt(alive*pz)]

  quat_rotate!(r, q)
  quat_rotate!(p, q)
  coords.q = vifelse(isnothing(coords.q), coords.q, quat_mul(q, coords.q))

  v[i,XI]  = vifelse(coords.state[i]==STATE_ALIVE, r[1], v[i,XI])
  v[i,YI]  = vifelse(coords.state[i]==STATE_ALIVE, r[2], v[i,YI])
  v[i,PXI] = vifelse(coords.state[i]==STATE_ALIVE, p[1], v[i,PXI])
  v[i,PYI] = vifelse(coords.state[i]==STATE_ALIVE, p[2], v[i,PYI])

  exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, ele_orient * L2 - r[3])
end

#

@makekernel fastgtpsa=true function track_alignment_bend!(i, coords::Coords, entering, 
          x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, tilde_m, beta_0, L)
  v = coords.v

end

