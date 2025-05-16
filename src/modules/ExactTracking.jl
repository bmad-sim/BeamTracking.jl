"""
    ExactTracking

Module implementing exact particle tracking through drifts and handling of misalignments.
"""

# Define the Exact tracking method
struct Exact end

# Number of temporaries needed for a single particle (number of columns in work matrix)
MAX_TEMPS(::Exact) = 1

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@makekernel function update_P0!(i, v, work, Brho_initial, Brho_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * Brho_initial / Brho_final
  end
  return v
end
=#

"""
    misalign!(i, v, work, x_offset, y_offset, sgn)

Apply misalignment offsets to particle coordinates.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `x_offset`: Horizontal offset
- `y_offset`: Vertical offset
- `sgn`: Sign (-1 for entering, 1 for exiting)
"""
# TODO: handle rotational misalignments
@makekernel function misalign!(i, v, work, x_offset, y_offset, sgn) #x_rot, y_rot, tilt,
  ##@assert sgn == -1 || sgn == 1 "Incorrect value for sgn (use -1 if entering, 1 if exiting)"
  @inbounds begin
    @FastGTPSA! v[i,XI] += sgn*x_offset
    @FastGTPSA! v[i,YI] += sgn*y_offset
  end
  return v
end

"""
    exact_drift!(i, v, work, L, tilde_m, gamsqr_0, beta_0)

Track a particle through a drift space using exact equations of motion.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `L`: Drift length
- `tilde_m`: Normalized mass
- `gamsqr_0`: Square of reference gamma
- `beta_0`: Reference beta
"""
@makekernel function exact_drift!(i, v, work, L, tilde_m, gamsqr_0, beta_0)
  ##@assert size(work, 2) >= 1 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 1) for exact_drift!"
  @inbounds begin @FastGTPSA! begin
    work[i,1] = sqrt((1.0 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
    v[i,XI]   = v[i,XI] + v[i,PXI] * L / work[i,1]
    v[i,YI]   = v[i,YI] + v[i,PYI] * L / work[i,1]
    v[i,ZI]   = v[i,ZI] - ( (1.0 + v[i,PZI]) * L
                  * ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                  / ( beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) * work[i,1]
                      * (beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) + work[i,1])
                    )
                )
  end end
  return v
end

end