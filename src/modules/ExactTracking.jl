#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Exact end

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
  #@assert size(work, 2) >= 1 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 1) for exact_drift!"
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

function sincu(x)
  return sinc(x/pi)
end

"""
    exact_bend!(i, v, work, theta, gtot, g, L, mc2, p0c, beta_ref) 

Tracks a particle through a sector bend via exact tracking. (no edge angles)

#Arguments
- 'theta'    -- 'g' * 'L'
- 'gtot'     -- 'g' + 'dg'
- 'p0c'      -- reference momentum in eV
- 'beta_ref' -- 'p0c' / sqrt('mc2'^2 + 'p0c'^2)
"""
@inline function exact_bend!(i, v, work, theta, gtot, g, L, mc2, p0c, beta_ref)
  @inbounds begin #@FastGTPSA! begin
      work[i,1] = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2) #pt
      work[i,2] = theta + asin(v[i,PXI] / work[i,1]) #phi1
      work[i,3] = gtot / work[i,1] #gp
      work[i,4] = 1+g*v[i,XI] # 1 + g * x
      work[i,5] = cos(work[i,2]) #cos(theta + phi1)
      work[i,6] = sincu(theta) #sincu(theta)
      work[i,7] = 2*work[i,4]*sin(work[i,2])*L*work[i,6]- work[i,3]*work[i,4]^2*L^2*(work[i,6])^2 #alpha
      if abs(work[i,2]) < π/2
          work[i,8] = work[i,7]/(sqrt(work[i,5]^2 + work[i,3]*work[i,7]) + work[i,5]) #xi
      else
          work[i,8] = (sqrt(work[i,5]^2 + work[i,3]*work[i,7]) - work[i,5]) / work[i,3] #xi
      end
      work[i,9] = -L*work[i,6]-v[i,XI]*sin(theta) #Lcv
      work[i,10] = 2 * (work[i,2] - atan(work[i,8], -work[i,9])) #theta_p
      work[i,11] = sqrt(work[i,9]^2 + work[i,8]^2) / sincu(work[i,10]/2) #Lp
      work[i,12] = (v[i,PZI] * p0c + p0c) #p

      v[i,XI] = v[i,XI]*cos(theta) - g/2*L^2*(sincu(theta/2))^2 + work[i,8]
      v[i,PXI] = work[i,1]*sin(work[i,2] - work[i,10])
      v[i,YI] = v[i,YI] + v[i,PYI]*work[i,11]/work[i,1] 
      v[i,ZI] = v[i,ZI] - (1 + v[i,PZI])*work[i,11]/work[i,1] + 
                      L*(work[i,12]/sqrt(mc2^2+work[i,12]^2))/beta_ref
  end #end
  return v
end

end