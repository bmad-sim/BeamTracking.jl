#---------------------------------------------------------------------------------------------------
# coord_concatenation(r1, q1, r2, q2) -> dr, q

"""
    coord_concatenation(r1, q1, r2, q2) -> dr, q

Returns the composite transform for the coordinate transform (`r1`, `q1`) followed by (`r2`, `q2`).
Note: This is the transform for coordinates which is the reverse of particle propagation.

## Output
- `r`     Composite coordinate origin shift.
- `q`     Composite Quaternion rotation.
"""
function coord_concatenation(r1, q1, r2, q2)
  return quat_rotate(r2, q1) .+ r1, quat_mul(q1, q2)
end



#---------------------------------------------------------------------------------------------------
# coord_alignment_bend_mid

"""
    coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, L)
                                                                  -> mid_r, mid_q, st, ct

Returns the `s`-independent part of the bend alignment transform: the misalignment
transform referenced to the bend arc center, plus the precomputed `sincos(tilt_ref)`.
Compute once, then feed to `coord_alignment_bend_at_s` for as many `s` values as you like.

## Returns
- `mid_r`   Arc-center-referenced origin shift.
- `mid_q`   Arc-center-referenced quaternion rotation.
- `st`      `sin(tilt_ref)`.
- `ct`      `cos(tilt_ref)`.
"""
@inline function coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, L)
  st, ct = sincos(tilt_ref)

  ang2 = L * g_ref / 2
  one_c = one_cos_norm(ang2)
  f = -L^2 * g_ref * one_c / 4

  # Translate from arc center to chord midpoint (q is still identity here).
  r = (f*ct, f*st, 0)

  # Misalignment transform (about the chord midpoint).
  r = r .+ (x_off, y_off, z_off)
  q = rot_quaternion(x_rot, y_rot, tilt)

  # Rotate by -tilt_ref.
  q = quat_mul(q, rot_quaternion(0, 0, tilt_ref))

  # Translate from chord midpoint back to arc center.
  r = r .+ quat_rotate((L^2 * g_ref * one_c / 4, 0, 0), q)

  return r, q, st, ct
end

#---------------------------------------------------------------------------------------------------
# coord_alignment_bend_at_s

"""
    coord_alignment_bend_at_s(mid_r, mid_q, st, ct, g_ref, L, s) -> r, q

Returns the bend alignment transform at arc-length `s` (0 ≤ s ≤ L) along the bend,
given the precomputed `(mid_r, mid_q, st, ct)` from `coord_alignment_bend_mid`.
At `s = 0` this is the entrance-face transform; at `s = L` it is the exit-face
transform (before inversion).

## Returns
- `r`     Coordinate origin shift.
- `q`     Quaternion rotation.
"""
@inline function coord_alignment_bend_at_s(mid_r, mid_q, st, ct, g_ref, L, s)
  @FastGTPSA begin
    ds  = L/2 - s
    ang = ds * g_ref
    rx  = -ds * ang * one_cos_norm(ang)
    rz  =  ds * sincu(ang)
  end
  sa, ca = sincos(ang/2)

  # First arc transform (ds, tilt_ref): tilt rotation reduces to a 2D x–y rotation.
  ra = (rx*ct, rx*st, rz)
  qa = (ca, sa*st, -sa*ct, 0)        # = rot_quat((st, -ct, 0), ang)

  # Second arc transform (-ds, 0): trq is identity ⇒ r unrotated; rz negates, rx shared.
  rb = (rx, 0, -rz)
  qb = (ca, 0, sa, 0)              # = rot_quat((0, -1, 0), -ang)

  # Concatenation 1: (ra, qa) ∘ (mid_r, mid_q)
  Q1 = quat_mul(qa, mid_q)
  R1 = quat_rotate(mid_r, qa) .+ ra

  # Concatenation 2: (R1, Q1) ∘ (rb, qb)
  q = quat_mul(Q1, qb)
  r = quat_rotate(rb, Q1) .+ R1
  return r, q
end

#---------------------------------------------------------------------------------------------------
#= coord_alignment_bend_entering

"""
    coord_alignment_bend_entering(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                              g_ref, tilt_ref, ele_orient, L) -> dr, q

Returns `dr` origin shift and `q` quaternion rotation for the coordinate transformation
from the nominal bend entrance face (in branch coordinates) to the actual entrance face
(in body coordinates) taking into account the element alignment parameters.
"""
@inline function coord_alignment_bend_entering(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                               g_ref, tilt_ref, ele_orient, L)
  mid_r, mid_q, st, ct = coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                                  g_ref, tilt_ref, L)
  return coord_alignment_bend_at_s(mid_r, mid_q, st, ct, g_ref, L, zero(L))
end

#---------------------------------------------------------------------------------------------------
# coord_alignment_bend_exiting

"""
    coord_alignment_bend_exiting(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                              g_ref, tilt_ref, ele_orient, L) -> dr, q

Returns `dr` origin shift and `q` quaternion rotation for the coordinate transformation
from the actual exit face (in body coordinates) to the nominal bend exit face
(in branch coordinates) taking into account the element alignment parameters.
"""
@inline function coord_alignment_bend_exiting(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                              g_ref, tilt_ref, ele_orient, L)
  mid_r, mid_q, st, ct = coord_alignment_bend_mid(x_off, y_off, z_off, x_rot, y_rot, tilt,
                                                  g_ref, tilt_ref, L)
  dr, q = coord_alignment_bend_at_s(mid_r, mid_q, st, ct, g_ref, L, L)
  q_inv = quat_inv(q)
  return .-quat_rotate(dr, q_inv), q_inv
end
=#