"""
    alignment_drift!(i, coords, L)

Drift used in the alignment algorithm. 
The difference between this and the standard drift is that the reference time does not change
as the particle drifts from beginning to end.

## Arguments
- `L`:       element length, in meters
"""
@inline function alignment_drift!(i, coords::Coords, L)
  @FastGTPSA begin @inbounds begin 
    v = coords.v

    rel_p = 1 + v[i,PZI]
    P_t2 = v[i,PXI]*v[i,PXI] + v[i,PYI]*v[i,PYI]
    P_s2 = rel_p*rel_p - P_t2

    good_momenta = (P_s2 > 0)
    alive = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)

    P_s2_1 = one(P_s2)
    P_s = sqrt(vifelse(good_momenta, P_s2, P_s2_1))

    new_x = v[i,XI] + v[i,PXI] * L / P_s
    new_y = v[i,YI] + v[i,PYI] * L / P_s
    new_z = v[i,ZI] - rel_p * L / P_s

    v[i,XI] = vifelse(alive, new_x, v[i,XI])
    v[i,YI] = vifelse(alive, new_y, v[i,YI])
    v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
  end end
end

