#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Exact end

MAX_TEMPS(::Exact) = 8

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
using ..BeamTracking: C_LIGHT
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

# Misalignments (TO-DO: rotational misalignments)
@makekernel function misalign!(i, v, work, x_offset, y_offset, sgn) #x_rot, y_rot, tilt,
  #@assert sgn == -1 || sgn == 1 "Incorrect value for sgn (use -1 if entering, 1 if exiting)"
  @inbounds begin
    @FastGTPSA! v[i,XI] += sgn*x_offset
    @FastGTPSA! v[i,YI] += sgn*y_offset
  end
  return v
end

@makekernel function exact_drift!(i, v, work, L, tilde_m, gamsqr_0, beta_0)
  @assert size(work, 2) >= 2 && size(work, 1) >= size(v,1) "Size of work matrix must be at least ($(size(v,1)), 2) for exact_drift!"
  @inbounds begin @FastGTPSA! begin
    work[i,1] = 1 + v[i,PZI]                                 # rel_p
    work[i,2] = sqrt(work[i,1]^2 - v[i,PXI]^2 - v[i,PYI]^2)  # ps
    v[i,XI]   += v[i,PXI] * L / work[i,2]
    v[i,YI]   += v[i,PYI] * L / work[i,2]
    v[i,ZI]   -=  work[i,1] * L *
                    ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) / 
                    ( beta_0 * sqrt(work[i,1]^2 + tilde_m^2) * work[i,2] * 
                      (beta_0 * sqrt(work[i,1]^2 + tilde_m^2) + work[i,2]) )
  end end
  return v
end

@makekernel function exact_solenoid!(i, v, work, L, ks, tilde_m, gamsqr_0, beta_0)
  @assert size(work, 2) >= 8 && size(work, 1) >= size(v,1) "Size of work matrix must be at least ($(size(v,1)), 8) for exact_solenoid!"
  @inbounds begin @FastGTPSA! begin
    # Recurring variables
    work[i,1] = 1 + v[i,PZI]                                 # rel_p
    work[i,2] = sqrt(work[i,1]^2 - (v[i,PXI] + v[i,YI] * ks / 2)^2 - (v[i,PYI] - v[i,XI] * ks / 2)^2) # pr
    work[i,3] = sin(ks * L / work[i,2])                      # S
    work[i,4] = 1 + cos(ks * L / work[i,2])                  # Cp
    work[i,5] = 1 - cos(ks * L / work[i,2])                  # Cm
    # Temporaries
    work[i,6] = v[i,XI]                                      # x_0
    work[i,7] = v[i,PXI]                                     # px_0
    work[i,8] = v[i,YI]                                      # y_0
    # Update
    v[i,ZI]  -= work[i,1] * L * 
                  ((v[i,PXI] + v[i,YI] * ks / 2)^2 + (v[i,PYI] - v[i,XI] * ks / 2)^2 - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) / 
                  ( beta_0 * sqrt(work[i,1]^2 + tilde_m^2) * work[i,2] * (beta_0 * sqrt(work[i,1]^2 + tilde_m^2) + work[i,2]) )
    v[i,XI] = work[i,4] * work[i,6] / 2 + work[i,3] * (work[i,7] / ks + work[i,8] / 2) + work[i,5] * v[i,PYI] / ks
    v[i,PXI] = work[i,3] * (v[i,PYI] / 2 - ks * work[i,6] / 4) + work[i,4] * work[i,7] / 2 - ks * work[i,5] * work[i,8] / 4
    v[i,YI] = work[i,3] * (v[i,PYI] / ks - work[i,6] / 2) + (work[i,4] * work[i,8] - work[i,5] * work[i,7]) / ks
    v[i,PYI]  = ks * work[i,5] * work[i,6] / 4 - work[i,3] * (work[i,7] / 2 + ks * work[i,8] / 4) + work[i,4] * v[i,PYI] / 2
  end end
  return v
end

@makekernel function patch!(i, v, work, L, tilde_m, gamsqr_0, beta_0, dt, dx, dy, dz, winv::Union{AbstractArray,Nothing})
  @assert size(work,2) >= 8 && size(work, 1) >= size(v, 1) "Size of work array must be at least ($(size(v, 1)), 9) for patch transformations. Received $work"
  @assert isnothing(winv) || (size(winv,1) == 3 && size(winv,2) == 3) "The inverse rotation matrix must be either `nothing` or 3x3 for patch!. Received $winv"
  @inbounds begin @FastGTPSA! begin
    # Temporary momentum [δp, pz]
    work[i,1] = 1 + v[i,PZI]                                 # rel_p
    work[i,2] = sqrt(work[i,1]^2 - v[i,PXI]^2 - v[i,PYI]^2)  # ps_0
  end end
    # Only apply rotations if needed
    if isnothing(winv)
      @inbounds begin @FastGTPSA! begin
      # No rotation case
      v[i,XI] -= dx
      v[i,YI] -= dy
      
      # Apply t_offset
      v[i,ZI] += work[i,1]*sqrt(work[i,1]/(work[i,1]^2+tilde_m^2))*C_LIGHT*dt

      # Drift to face
      v = exact_drift!(i, v, work, -dz, tilde_m, gamsqr_0, beta_0)
      end end
    else
      @inbounds begin @FastGTPSA! begin
      # Translate position vector [x, y]
      work[i,3] = v[i,XI] - dx                                # x_0
      work[i,4] = v[i,YI] - dy                                # y_0

      # Temporary momentum vector [px, py]
      work[i,5] = v[i,PXI]                                    # px_0
      work[i,6] = v[i,PYI]                                    # py_0
      
      # Transform position vector [x - dx, y - dy, -dz]
      v[i,XI]   = winv[1,1]*work[i,3] + winv[1,2]*work[i,4] - winv[1,3]*dz
      v[i,YI]   = winv[2,1]*work[i,3] + winv[2,2]*work[i,4] - winv[2,3]*dz
      work[i,7] = winv[3,1]*work[i,3] + winv[3,2]*work[i,4] - winv[3,3]*dz  # s_f
      
      # Transform momentum vector [px, py, ps]
      v[i,PXI]  = winv[1,1]*work[i,5] + winv[1,2]*work[i,6] + winv[1,3]*work[i,2]
      v[i,PYI]  = winv[2,1]*work[i,5] + winv[2,2]*work[i,6] + winv[2,3]*work[i,2]
      work[i,8] = winv[3,1]*work[i,5] + winv[3,2]*work[i,6] + winv[3,3]*work[i,2] # ps_f

      # Apply t_offset
      v[i,ZI] += work[i,1]/sqrt(work[i,1]^2+tilde_m^2)*C_LIGHT*dt

      # Drift to face
      v[i,XI] -= work[i,7] * v[i,PXI] / work[i,8]
      v[i,YI] -= work[i,7] * v[i,PYI] / work[i,8]
      v[i,ZI] += work[i,7] * work[i,1] / work[i,8] + L*work[i,1]*sqrt((1+tilde_m^2)/(work[i,1]^2+tilde_m^2))
      end end
    end
  return v
end


# Utility functions ============================================================

# Rotation matrix
"""
  w_matrix(x_rot, y_rot, z_rot)

Constructs a rotation matrix based on the given Bryan-Tait angles.

Bmad/SciBmad follows the MAD convention of applying z, x, y rotations in that order.
Furthermore, in ReferenceFrameRotations, the rotation angles around Y and Z axes 
are defined as negative of the SciBmad `y_rot` and `z_rot`.

The inverse matrix reverses the order of operations and their signs.


Arguments:
- `x_rot::Real`: Rotation angle around the x-axis.
- `y_rot::Real`: Rotation angle around the y-axis.
- `z_rot::Real`: Rotation angle around the z-axis.

Returns:
- `Matrix{Float64}`: Rotation matrix.
"""
function w_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(-z_rot, x_rot, -y_rot, :ZXY)
end

# Inverse rotation matrix
function w_inv_matrix(x_rot, y_rot, z_rot)
  return ReferenceFrameRotations.angle_to_rot(y_rot, -x_rot, z_rot, :YXZ)
end
end