#---------------------------------------------------------------------------------------------------
# coord_concatenation(r1, q1, r2, q2) -> dr, q

"""
    coord_concatenation(r1, q1, r2, q2) -> dr, q

Returns the composite transform for the coordinate transform (`r1`, `q1`) followed by (`r2`, `q2`).

## Output
- `r`     Composite coordinate origin shift.
- `q`     Composite Quaternion rotation.
"""
function coord_concatenation(r1, q1, r2, q2)
  return quat_rotate(r2, q1) + r1, quat_mul(q1, q2)
end

#---------------------------------------------------------------------------------------------------
# coord_bend_arc_transform

"""
    function coord_bend_arc_transform(ds, g_ref, tilt_ref) -> r, q

Return the transform to the coordinate system that is rotated about the bend axis by an
angle `ds*g_ref`.

## Arguments
- `ds`        Distance along the bend arc the coordinate system is transformed.
- `g_ref`     Reference `g = 1/radius`.
- `tilt_ref`  Reference tilt.

## Output
- `r`     Coordinate origin shift.
- `q`     Quaternion rotation.
"""
@inline function coord_bend_arc_transform(ds, g_ref, tilt_ref)
  @FastGTPSA begin @inbounds begin 
    ang = ds * g_ref
    axis = (sin(tilt_ref), -cos(tilt_ref), 0.0)
    q = rot_quat(axis, ang)

    trq = rot_quat((0.0, 0.0, 1.0), tilt_ref)
    r = (-ds * ang * one_cos_norm(ang), 0.0, ds * sincu(ang))
    r = quat_rotate(r, trq)

    return r, q
  end end
end

#---------------------------------------------------------------------------------------------------

"""
    coord_alignment_bend_entering(x_off, y_off, z_off, x_rot, y_rot, tilt, 
                                              g_ref, tilt_ref, ele_orient, L) -> dr, q

Returns `dr` origin shift and `q` quaternion rotation for the coordinate transformation
from the nominal bend entrance face (in branch coordinates) to the actual entrance face
(in body coordinates) taking into account the element alignment parameters.

## Arguments
- `x_off`, `y_off`, `z_off`   Element offset.
- `x_rot`, `y_rot`, 'tilt`    Element orientation.
- `g_ref`                     Reference g = 1/bend radius.
- `tilt_ref`                  Branch coords tilt.
- `ele_orient`                Element longitudinal orientation: +1 => normal, -1 => reversed.
- `L`                         Element length.

## Returns
- `dr`    Coordinate origin shift.
- `q`     Quaternion rotation.
"""

@inline function coord_alignment_bend_entering(x_off, y_off, z_off, 
                                        x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)

  # Start: Transform to coords at center of bend arc.
  r, q = coord_bend_arc_transform(0.5*L, g_ref, tilt_ref)
  ## println("***A: $r  :: $q")

  # Translate to coords with origin at element chord midpoint.
  ang2 = 0.5 * L * g_ref
  f = -0.25 * L^2 * g_ref * one_cos_norm(ang2)
  r += quat_rotate((f*cos(tilt_ref), f*sin(tilt_ref), 0.0), q)
  ## println("***B: $r  :: $q")

  # Misalignment transform
  r += quat_rotate((x_off, y_off, z_off), q) 
 
  dq = rot_quaternion(x_rot, y_rot, tilt)
  q = quat_mul(q, dq)
  ## println("***W: $r  :: $q")

  # Rotating by -tilt_ref
  dq = rot_quaternion(0.0, 0.0, tilt_ref)
  q = quat_mul(q, dq)
  ## println("***X: $r  :: $q")

  # Translate from chord center to arc center.
  dr = (0.25 * L^2 * g_ref * one_cos_norm(ang2), 0.0, 0.0)
  r += quat_rotate(dr, q)
  ## println("***Y: $r  :: $q")

  # Transform from arc center back to entrance face.
  dr, dq = coord_bend_arc_transform(-0.5*L, g_ref, 0.0)
  r, q = coord_concatenation(r, q, dr, dq)
  ## println("***Z: $r  :: $q")
  return r, q
end
 
#---------------------------------------------------------------------------------------------------

@inline function coord_alignment_bend_exiting(x_off, y_off, z_off, 
                                        x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)

  # Idea: Transform from non-misaligned (branch) coords to body coords and then
  # return the inverse.

  # Start: Transform to coords at center of bend arc.
  r, q = coord_bend_arc_transform(-0.5*L, g_ref, tilt_ref)
  ## println("***A: $r  :: $q")

  # Translate to coords with origin at element chord midpoint.
  ang2 = 0.5 * L * g_ref
  f = -0.25 * L^2 * g_ref * one_cos_norm(ang2)
  r += quat_rotate((f*cos(tilt_ref), f*sin(tilt_ref), 0.0), q)
  ## println("***B: $r  :: $q")

  # Misalignment transform
  r += quat_rotate((x_off, y_off, z_off), q) 
 
  dq = rot_quaternion(x_rot, y_rot, tilt)
  q = quat_mul(q, dq)
  ## println("***W: $r  :: $q")

  # Rotating by -tilt_ref
  dq = rot_quaternion(0.0, 0.0, tilt_ref)
  q = quat_mul(q, dq)
  ## println("***X: $r  :: $q")

  # Translate from chord center to arc center.
  dr = (0.25 * L^2 * g_ref * one_cos_norm(ang2), 0.0, 0.0)
  r += quat_rotate(dr, q)
  ## println("***Y: $r  :: $q")

  # Transform from arc center back to entrance face.
  dr, dq = coord_bend_arc_transform(0.5*L, g_ref, 0.0)
  r, q = coord_concatenation(r, q, dr, dq)
  ## println("***Z: $r  :: $q")

  # Return inverse transform
  q = quat_inv(q)
  return -quat_rotate(r, q), q
end
