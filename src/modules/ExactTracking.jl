#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Exact end

MAX_TEMPS(::Exact) = 1

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@inline function update_P0!(i, v, work, Brho_initial, Brho_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * Brho_initial / Brho_final
  end
  return v
end
=#

# Misalignments (TO-DO: rotational misalignments)
@inline function misalign!(i, v, work, x_offset, y_offset, sgn) #x_rot, y_rot, tilt,
  @assert sgn == -1 || sgn == 1 "Incorrect value for sgn (use -1 if entering, 1 if exiting)"
  @inbounds begin
    @FastGTPSA! v[i,XI] += sgn*x_offset
    @FastGTPSA! v[i,YI] += sgn*y_offset
  end
  return v
end

@inline function exact_drift!(i, v, work, L, tilde_m, gamsqr_0, beta_0)
  @assert size(work, 2) >= 1 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 1) for exact_drift!"
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

@inline function patch!(i, v, work, dx, dy, dz, y_rot, x_rot, z_rot)
  @inbounds begin @FastGTPSA! begin
    work[i,1] = 1 + v[i,PZI]  # rel_p
    work[i,2] = sqrt(work[i,1]^2 - v[i,PXI]^2 - v[i,PYI]^2)  # pz
    
    # Apply rotations only if needed
    if x_rot != 0 || y_rot != 0 || z_rot != 0
      work[i,3] = cos(z_rot)
      work[i,4] = sin(z_rot)
      work[i,5] = cos(x_rot)
      work[i,6] = sin(x_rot)
      work[i,7] = cos(y_rot)
      work[i,8] = sin(y_rot)
      
      work[i,9] = v[i,XI] - dx
      work[i,10] = v[i,YI] - dy
      
      # Transform position vector [x - dx, y - dy, -dz]
      work[i,11] = (work[i,7]*work[i,3])*work[i,9] + 
                   (work[i,7]*work[i,4])*work[i,10] + 
                   (-work[i,8])*(-dz)
                   
      work[i,12] = (work[i,6]*work[i,8]*work[i,3] - work[i,5]*work[i,4])*work[i,9] + 
                   (work[i,6]*work[i,8]*work[i,4] + work[i,5]*work[i,3])*work[i,10] + 
                   (work[i,6]*work[i,7])*(-dz)
                   
      work[i,13] = (work[i,5]*work[i,8]*work[i,3] + work[i,6]*work[i,4])*work[i,9] + 
                   (work[i,5]*work[i,8]*work[i,4] - work[i,6]*work[i,3])*work[i,10] + 
                   (work[i,5]*work[i,7])*(-dz)
      
      work[i,14] = v[i,PXI]
      work[i,15] = v[i,PYI]
      
      # Transform momentum vector [px, py, pz]
      v[i,PXI] = (work[i,7]*work[i,3])*work[i,14] + 
                 (work[i,7]*work[i,4])*work[i,15] + 
                 (-work[i,8])*work[i,2]
                 
      v[i,PYI] = (work[i,6]*work[i,8]*work[i,3] - work[i,5]*work[i,4])*work[i,14] + 
                 (work[i,6]*work[i,8]*work[i,4] + work[i,5]*work[i,3])*work[i,15] + 
                 (work[i,6]*work[i,7])*work[i,2]
                 
      work[i,16] = (work[i,5]*work[i,8]*work[i,3] + work[i,6]*work[i,4])*work[i,14] + 
                   (work[i,5]*work[i,8]*work[i,4] - work[i,6]*work[i,3])*work[i,15] + 
                   (work[i,5]*work[i,7])*work[i,2]
      
      v[i,XI] = work[i,11]
      v[i,YI] = work[i,12]
      
      work[i,17] = (work[i,5]*work[i,8]*work[i,3] + work[i,6]*work[i,4])*dx + 
                   (work[i,5]*work[i,8]*work[i,4] - work[i,6]*work[i,3])*dy + 
                   (work[i,5]*work[i,7])*dz

      # Drift to face
      v[i,XI] -= work[i,13] * v[i,PXI] / work[i,16]
      v[i,YI] -= work[i,13] * v[i,PYI] / work[i,16]
      v[i,ZI] += work[i,13] * work[i,1] / work[i,16] + work[i,17]
    else
      # No rotation case
      v[i,XI] -= dx
      v[i,YI] -= dy
      
      # Drift to face
      v[i,XI] += dz * v[i,PXI] / work[i,2]
      v[i,YI] += dz * v[i,PYI] / work[i,2]
      v[i,ZI] -= dz * work[i,1] / work[i,2] + work[i,17]
    end
  end end
  
  return v
end
end