@makekernel fastgtpsa=true function track_alignment_straight_entering!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  q = inv_rot_quaternion(x_rot, y_rot, tilt)
  L2 = 0.5 * L * ele_orient

  translation!(i, coords, (x_off, y_off, 0.0), 0.0)
  dz_new = rotation!(i, coords, q, -L2 - z_off)
  isochronous_drift!(i, coords, -dz_new - L2)
end

#

#---------------------------------------------------------------------------------------------------

@makekernel fastgtpsa=true function track_alignment_straight_exiting!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  v = coords.v
  L2 = 0.5 * L * ele_orient

  q = rot_quaternion(x_rot, y_rot, tilt)
  dz_new = rotation!(i, coords, q, L2)

  alive = (coords.state[i] == STATE_ALIVE)
  v[i,XI] = vifelse(alive, v[i,XI] + x_off, v[i,XI])
  v[i,YI] = vifelse(alive, v[i,YI] + y_off, v[i,YI])

  isochronous_drift!(i, coords, L2 - z_off - dz_new)
end

#---------------------------------------------------------------------------------------------------
"""
    function track_coord_transform!(i, coords::Coords, r, q)

Particle tracking using the coordinate transformation `(r, q)` followed by a drift to
the "element edge". The drifting is such that the particle's `z` coordinate (not to be confused with
the phase-space `z` coordinate) is zero.

`(r, q)` is the transformation of the coordinate system. The particle transformation is the
inverse of this transformation.
"""

@makekernel fastgtpsa=true function track_coord_transform!(i, coords::Coords, r, q)
  translation!(i, coords, (r[1], r[2], 0.0), 0.0)
  dz_new = rotation!(i, coords, quat_inv(q), -r[3])
  ##println("***AA: $(coords.v[1,:]) :: $dz_new")
  isochronous_drift!(i, coords, -dz_new)
  ##println("***BB: $(coords.v[1,:])")
end

