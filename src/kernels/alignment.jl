#---------------------------------------------------------------------------------------------------
# track_alignment_straight_at_s!
# in=true  : nominal -> body  (entering),  reduces to track_alignment_straight_entering! at s=0
# in=false : body -> nominal  (exiting),   reduces to track_alignment_straight_exiting!  at s=L

@inline function track_alignment_straight_at_s!(i, coords::Coords, x_off, y_off, z_off,
                                  x_rot, y_rot, tilt, ele_orient, L, s, ::Val{in}) where {in}
  # Signed longitudinal distance from the element center (rotation center) to the s-face.
  L2 = @FastGTPSA (0.5*L - s) * ele_orient

  if in
    q = inv_rot_quaternion(x_rot, y_rot, tilt)
    translation!(i, coords, (x_off, y_off, 0), 0)
    dz_new = rotation!(i, coords, q, -L2 - z_off)
    isochronous_drift!(i, coords, -dz_new - L2)
  else
    q = rot_quaternion(x_rot, y_rot, tilt)
    dz_new = rotation!(i, coords, q, -L2)
    translation!(i, coords, (-x_off, -y_off, 0), 0)     # adds x_off, y_off (guarded)
    isochronous_drift!(i, coords, -L2 - z_off - dz_new)
  end
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
@makekernel function track_coord_transform!(i, coords::Coords, r, q)
  translation!(i, coords, (r[1], r[2], 0), 0)
  dz_new = rotation!(i, coords, quat_inv(q), -r[3])
  ##println("***AA: $(coords.v[1,:]) :: $dz_new")
  isochronous_drift!(i, coords, -dz_new)
  ##println("***BB: $(coords.v[1,:])")
end

#---------------------------------------------------------------------------------------------------
# s = 0 entering, s = L exiting
@inline function track_coord_bend_transform_at_s!(i, coords::Coords, mid_r, mid_q, st, ct, g_ref, L, s, ::Val{in}) where {in}
  r, q = coord_alignment_bend_at_s(mid_r, mid_q, st, ct, g_ref, L, s)
  if in
    track_coord_transform!(i, coords, r, q)
  else
    q_inv = quat_inv(q)
    r_inv = .-quat_rotate(r, q_inv)
    track_coord_transform!(i, coords, r_inv, q_inv)
  end
end