"""
    function coord_bend_transform!(i, coords::Coords, ds, g_ref, tilt_ref, z0) -> z1

Transform phase-space particle coordinates when the coordinate system is transformed along
a bend arc a distance `ds`.

## Arguments
- `ds`:       Distance along the bend arc the coordinate system is transformed.
- `g_ref`     Reference `g = 1/radius`.
- `tilt_ref`  Reference tilt.
-  `z0`       The particle longitudinal distance from the center of rotation. 

## Output
- `z1`        The particle longitudinal distance from the center of rotation after rotation. 

"""
@inline function coord_bend_transform!(i, coords::Coords, ds, g_ref, tilt_ref, z0)
  @FastGTPSA begin @inbounds begin 
    v = coords.v

    angle = ds * g_ref
    L = ds * [-0.5 * angle * sincu(0.5*angle), 0.0, sincu(angle)]
    q = [0.0]

  end end
end