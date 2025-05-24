#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Exact end

MAX_TEMPS(::Exact) = 9

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
    work[i,1] = 1 + v[i,PZI]                                 # 1+δp
    work[i,2] = sqrt(work[i,1]^2 - v[i,PXI]^2 - v[i,PYI]^2)  # P_s
    v[i,XI]   += v[i,PXI] * L / work[i,2]
    v[i,YI]   += v[i,PYI] * L / work[i,2]
    v[i,ZI]   -=  work[i,1] * L *
                    ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0) / 
                    ( beta_0 * sqrt(work[i,1]^2 + tilde_m^2) * work[i,2] * 
                      (beta_0 * sqrt(work[i,1]^2 + tilde_m^2) + work[i,2])
                    )
  end end
  return v
end

@makekernel function exact_solenoid!(i, v, work, L, ks, p0c, mc2)
  @assert size(work, 2) >= 8 && size(work, 1) >= size(v,1) "Size of work matrix must be at least ($(size(v,1)), 8) for exact_solenoid!"
  @inbounds begin @FastGTPSA! begin
    # Recurring variables
    work[i,1] = 1 + v[i,PZI]
    work[i,2] = sqrt(work[i,1]^2 - (v[i,PXI] + v[i,YI] * ks / 2)^2 - (v[i,PYI] - v[i,XI] * ks / 2)^2)
    work[i,3] = sin(ks * L / work[i,2])
    work[i,4] = 1 + cos(ks * L / work[i,2])
    work[i,5] = 1 - cos(ks * L / work[i,2])
    # Temporaries
    work[i,6] = work[i,4] * v[i,XI] / 2 + work[i,3] * (v[i,PXI] / ks + v[i,YI] / 2) + work[i,5] * v[i,PYI] / ks
    work[i,7] = work[i,3] * (v[i,PYI] / 2 - ks * v[i,XI] / 4) + work[i,4] * v[i,PXI] / 2 - ks * work[i,5] * v[i,YI] / 4
    work[i,8] = work[i,3] * (v[i,PYI] / ks - v[i,XI] / 2) + (work[i,4] * v[i,YI] - work[i,5] * v[i,PXI]) / ks
    # Update
    v[i,PYI]  = ks * work[i,5] * v[i,XI] / 4 - work[i,3] * (v[i,PXI] / 2 + ks * v[i,YI] / 4) + work[i,4] * v[i,PYI] / 2
    v[i,XI]   = work[i,6]
    v[i,PXI]  = work[i,7]
    v[i,YI]   = work[i,8]
    v[i,ZI]  += work[i,1] * L * ( sqrt(p0c^2 + mc2^2) / sqrt((p0c*work[i,1])^2 + mc2^2) - 1 / work[i,2] )
  end end
  return v
end

@makekernel function patch!(i, v, work, L, p0c, mc2, dt, dx, dy, dz, winv::Union{AbstractArray,Nothing})
  @assert size(work,2) >= 9 && size(work, 1) >= size(v, 1) "Size of work array must be at least ($(size(v, 1)), 9) for patch transformations. Received $work"
  @assert isnothing(winv) || (size(winv,1) == 3 && size(winv,2) == 3) "The inverse rotation matrix must be either `nothing` or 3x3 for patch!. Received $winv"
  @inbounds begin @FastGTPSA! begin
    # Temporary momentum [δp, pz]
    work[i,1] = 1 + v[i,PZI]                                 # 1+δp
    work[i,2] = sqrt(work[i,1]^2 - v[i,PXI]^2 - v[i,PYI]^2)  # p_s
  end end
    # Only apply rotations if needed
    if isnothing(winv)
      @inbounds begin @FastGTPSA! begin
      # No rotation case
      v[i,XI] -= dx
      v[i,YI] -= dy
      
      # Apply t_offset
      v[i,ZI] += work[i,1]*sqrt(p0c*work[i,1]/((p0c*work[i,1])^2+mc2^2))*C_LIGHT*dt

      # Drift to face
      v = exact_drift!(i, v, work, -dz, p0c, mc2)
      end end
    else
      @inbounds begin @FastGTPSA! begin
      # Translate position vector [x, y]
      work[i,3] = v[i,XI] - dx
      work[i,4] = v[i,YI] - dy

      # Temporary momentum vector [px, py]
      work[i,5] = v[i,PXI]
      work[i,6] = v[i,PYI]
      
      # Transform position vector [x - dx, y - dy, -dz]
      v[i,XI]   = winv[1,1]*work[i,3] + winv[1,2]*work[i,4] - winv[1,3]*dz
      v[i,YI]   = winv[2,1]*work[i,3] + winv[2,2]*work[i,4] - winv[2,3]*dz
      work[i,7] = winv[3,1]*work[i,3] + winv[3,2]*work[i,4] - winv[3,3]*dz
      
      # Transform momentum vector [px, py, pz]
      v[i,PXI]  = winv[1,1]*work[i,5] + winv[1,2]*work[i,6] + winv[1,3]*work[i,2]
      v[i,PYI]  = winv[2,1]*work[i,5] + winv[2,2]*work[i,6] + winv[2,3]*work[i,2]
      work[i,8] = winv[3,1]*work[i,5] + winv[3,2]*work[i,6] + winv[3,3]*work[i,2]

      # Apply t_offset
      v[i,ZI] += p0c*work[i,1]/sqrt((p0c*work[i,1])^2+mc2^2)*C_LIGHT*dt

      # Drift length
      work[i,9] = winv[3,1]*dx + winv[3,2]*dy + winv[3,3]*dz

      # Drift to face
      v[i,XI] -= work[i,7] * v[i,PXI] / work[i,8]
      v[i,YI] -= work[i,7] * v[i,PYI] / work[i,8]
      v[i,ZI] += work[i,7] * work[i,1] / work[i,8] + L*work[i,1]*sqrt((p0c^2+mc2^2)/((p0c*work[i,1])^2+mc2^2))
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