# BeamTracking

Documentation currently in development. See [`SciBmad.jl`](https://github.com/bmad-sim/SciBmad.jl) for examples of how to track in the `SciBmad` accelerator physics software ecosystem.


## Reference time for an element that has a finite misalignment. 

By convention, when a particle's position is expressed in element body coordinates,
the reference time at any point in an element is independent of the element alignment. 
This convention simplifies the conversion between phase space `z` and the absolute time since
the conversion is independent of the alignment and only dependent upon where the particle is
with respect to the element being tracked through.

### Alignment transformation and 

The reference time convention discussed above has the following consequence: 
In the transformation from a particle at the 
nominal (no alignment) element edge, with the particle's position expressed in branch coordinates,
to the particle at the the actual element edge with the particle's position expressed in
body coordinates, the particle must be "drifted" from the starting position to the actual element face.
However, the reference time at the starting position, with branch coordinates, is the
same as the ending position with body coordinates. Thus the drifting transformation
of phase space `z` is different than the standard `z` transformation for an element that is a drift. 
A similar situation happens when the particle exits the element and there is a transformation
from body coordinate to branch coordinates along with a "drift" to the nominal downstream edge.


## RF phase and absolute time

By default, RF cavities use relative-time tracking: `phi0` is the phase seen by the
reference particle at each cavity, independent of the bunch reference time. Pass
`absolute_time_tracking=true` to `track!` to phase-lock every cavity to the global
reference time carried by `bunch.t_ref`. In this mode the phase seen by a particle is
shifted by `2*pi*rf_frequency*bunch.t_ref`, enabling effects caused by changes in the
physical RF frequency, such as radial steering.

The default remains relative-time tracking for compatibility with lattices whose cavity
phases are specified locally.