@inline function alignment(tm, bunch, alignmentparams, bendparams, L, entering::Bool)
  if !isactive(alignmentparams); return nothing; end

  x_off = alignmentparams.x_offset
  y_off = alignmentparams.y_offset
  z_off = alignmentparams.z_offset
  x_rot = alignmentparams.x_rot
  y_rot = alignmentparams.y_rot
  tilt  = alignmentparams.tilt

  ele_orient = 1   ## Future work: Need to extend this for reversed elements.

  #

  if isactive(bendparams) && (bendparams.g_ref != 0 || bendparams.tilt_ref != 0)
    if entering
      dr, q = coord_transform_alignment_bend_entering(x_off, y_off, z_off, 
                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, L)
      return KernelCall(BeamTracking.track_coord_transform!, (dr, q))
    else
      dr, q = coord_transform_alignment_bend_exiting(x_off, y_off, z_off, 
                x_rot, y_rot, tilt, bendparams.g_ref, bendparams.tilt_ref, ele_orient, L)
      return KernelCall(BeamTracking.track_coord_transform!, (dr, q))
    end

  else
    if entering
      return KernelCall(BeamTracking.track_alignment_straight_entering!, (x_off, y_off, z_off, 
                                                     x_rot, y_rot, tilt, ele_orient, L))
    else
      return KernelCall(BeamTracking.track_alignment_straight_exiting!, (x_off, y_off, z_off, 
                                                     x_rot, y_rot, tilt, ele_orient, L))
    end
  end
end

#---------------------------------------------------------------------------------------------------
# coord_transform_alignment_bend_entering

"""
    coord_transform_alignment_bend_entering(x_off, y_off, z_off, x_rot, y_rot, tilt, 
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

@inline function coord_transform_alignment_bend_entering(x_off, y_off, z_off, x_rot, y_rot, tilt, 
                                                                     g_ref, tilt_ref, ele_orient, L)

    #(dr, q) = coord_bend_transform!(0.5*L, g_ref, tilt_ref)
    
end