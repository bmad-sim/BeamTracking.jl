"""
    LinearTracking

Module implementing linear particle tracking methods.

This module provides functions for linear particle tracking through elements,
including drifts, quadrupoles, solenoids, and bends, using first-order approximations.
"""

# Define the Linear tracking method
struct Linear end
# Number of temporaries needed for a single particle (number of columns in work matrix)
MAX_TEMPS(::Linear) = 5

module LinearTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
const TRACKING_METHOD = Linear

"""
    linear_drift!(i, v, work, L, r56)

Track a particle through a drift space using linear approximation.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `L`: Drift length
- `r56`: Longitudinal dispersion
"""
@makekernel function linear_drift!(i, v, work, L, r56)
  @inbounds begin @FastGTPSA! begin
    v[i,XI] += v[i,PXI] * L
    v[i,YI] += v[i,PYI] * L
    v[i,ZI] += v[i,PZI] * r56
  end end
  return v
end


#=
 Generic function for an uncoupled matrix with coasting plane:

[ mx      0       0   d[1:2]]
[ 0       my      0   d[3:4]]
[ t[1:2]  t[3:4]  1   r56   ]

=#
"""
    linear_coast_uncoupled!(i, v, work, mx, my, r56, d, t)

Track a particle through an uncoupled element

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `mx`: 2x2 horizontal transfer matrix
- `my`: 2x2 vertical transfer matrix
- `r56`: Momentum compaction term
- `d`: Dispersion vector (optional)
- `t`: Path length terms (optional)
"""
@makekernel function linear_coast_uncoupled!(i, v, work, mx::AbstractMatrix, my::AbstractMatrix, r56, d::Union{AbstractArray,Nothing}, t::Union{AbstractArray,Nothing})
  ##@assert size(work, 2) >= 1 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 1) for linear_coast_uncoupled!"
  ##@assert size(mx) == (2,2) "Size of matrix mx must be (2,2) for linear_coast_uncoupled!. Received $(size(mx))"
  ##@assert size(my) == (2,2) "Size of matrix my must be (2,2) for linear_coast_uncoupled!. Received $(size(my))"
  ##@assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast_uncoupled!. Received $d"
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
    end end
  end
  @inbounds begin @FastGTPSA! begin
    work[i,1]= v[i,XI]
    v[i,XI]  = mx[1,1] * v[i,XI]   + mx[1,2] * v[i,PXI] 
    v[i,PXI] = mx[2,1] * work[i,1] + mx[2,2] * v[i,PXI]
    work[i,1]= v[i,YI]
    v[i,YI]  = my[1,1] * v[i,YI]   + my[1,2] * v[i,PYI] 
    v[i,PYI] = my[2,1] * work[i,1] + my[2,2] * v[i,PYI]
    v[i,ZI] += r56 * v[i,PZI]
  end end
  if !isnothing(d)
    @inbounds begin @FastGTPSA! begin
      v[i,XI]  += d[XI]  * v[i,PZI]
      v[i,PXI] += d[PXI] * v[i,PZI]
      v[i,YI]  += d[YI]  * v[i,PZI]
      v[i,PYI] += d[PYI] * v[i,PZI]
    end end
  end
  return v
end

"""
    linear_coast!(i, v, work, mxy, r56, d, t)

Track a particle through a coupled element

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `mxy`: 4x4 coupled transfer matrix
- `r56`: Momentum compaction term
- `d`: Dispersion vector (optional)
- `t`: Path length terms (optional)
"""
@makekernel function linear_coast!(i, v, work, mxy::AbstractMatrix, r56, d::Union{AbstractArray,Nothing}, t::Union{AbstractArray,Nothing})
  #@assert size(work, 2) >= 3 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 3) for linear_coast!"
  #@assert size(mxy) == (4,4) "Size of matrix mxy must be (4,4) for linear_coast!. Received $(size(mxy))"
  #@assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast!. Received $d"
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
    end end
  end
  @inbounds begin @FastGTPSA! begin
    work[i,1]= v[i,XI]
    work[i,2]= v[i,PXI]
    work[i,3]= v[i,YI]
    v[i,XI]  = mxy[XI, XI] * v[i,XI]   + mxy[XI, PXI] * v[i,PXI]  + mxy[XI, YI] * v[i,YI]   + mxy[XI, PYI] * v[i,PYI]
    v[i,PXI] = mxy[PXI,XI] * work[i,1] + mxy[PXI,PXI] * v[i,PXI]  + mxy[PXI,YI] * v[i,YI]   + mxy[PXI,PYI] * v[i,PYI]
    v[i,YI]  = mxy[YI, XI] * work[i,1] + mxy[YI, PXI] * work[i,2] + mxy[YI, YI] * v[i,YI]   + mxy[YI, PYI] * v[i,PYI] 
    v[i,PYI] = mxy[PYI,XI] * work[i,1] + mxy[PYI,PXI] * work[i,2] + mxy[PYI,YI] * work[i,3] + mxy[PYI,PYI] * v[i,PYI]
    v[i,ZI] += r56 * v[i,PZI]
  end end
  if !isnothing(d)
    @inbounds begin @FastGTPSA! begin
      v[i,XI]  += d[XI]  * v[i,PZI]
      v[i,PXI] += d[PXI] * v[i,PZI]
      v[i,YI]  += d[YI]  * v[i,PZI]
      v[i,PYI] += d[PYI] * v[i,PZI]
    end end
  end
  return v
end

"""
    linear_6D!(i, v, work, m)

Track a particle using a full 6D transfer matrix.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix (must be at least size (N_particle, 5))
- `m`: 6x6 transfer matrix

# Notes
- Handles full 6D coupled motion
- Uses work matrix for temporary calculations
"""
@makekernel function linear_6D!(i, v, work, m::AbstractMatrix)
  ##@assert size(work, 2) >= 5 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 5) for linear_6D!"
  ##@assert size(m) == (6,6) "Size of matrix m must be (6,6) for linear_6D!. Received $(size(m))"
  @inbounds begin @FastGTPSA! begin
    work[i,1]= v[i,XI]
    work[i,2]= v[i,PXI]
    work[i,3]= v[i,YI]
    work[i,4]= v[i,PYI]
    work[i,5]= v[i,ZI]
    v[i,XI]  = m[XI, XI] * v[i,XI]   + m[XI, PXI] * v[i,PXI]  + m[XI, YI] * v[i,YI]   + m[XI, PYI] * v[i,PYI]  + m[XI, ZI] * v[i,ZI]   + m[XI, PZI] * v[i,PZI]
    v[i,PXI] = m[PXI,XI] * work[i,1] + m[PXI,PXI] * v[i,PXI]  + m[PXI,YI] * v[i,YI]   + m[PXI,PYI] * v[i,PYI]  + m[PXI,ZI] * v[i,ZI]   + m[PXI,PZI] * v[i,PZI]
    v[i,YI]  = m[YI, XI] * work[i,1] + m[YI, PXI] * work[i,2] + m[YI, YI] * v[i,YI]   + m[YI, PYI] * v[i,PYI]  + m[YI, ZI] * v[i,ZI]   + m[YI, PZI] * v[i,PZI]
    v[i,PYI] = m[PYI,XI] * work[i,1] + m[PYI,PXI] * work[i,2] + m[PYI,YI] * work[i,3] + m[PYI,PYI] * v[i,PYI]  + m[PYI,ZI] * v[i,ZI]   + m[PYI,PZI] * v[i,PZI]
    v[i,ZI]  = m[ZI, XI] * work[i,1] + m[ZI, PXI] * work[i,2] + m[ZI, YI] * work[i,3] + m[ZI, PYI] * work[i,4] + m[ZI, ZI] * v[i,ZI]   + m[ZI, PZI] * v[i,PZI]
    v[i,PZI] = m[PZI,XI] * work[i,1] + m[PZI,PXI] * work[i,2] + m[PZI,YI] * work[i,3] + m[PZI,PYI] * work[i,4] + m[PZI,ZI] * work[i,5] + m[PZI,PZI] * v[i,PZI]
  end end
end

# Utility functions to create a linear matrix
"""
    linear_quad_matrices(K1, L)

Generate transfer matrices for a thick quadrupole.

# Arguments
- `K1`: Quadrupole strength, focusing and defocusing matrices based on K1 sign
- `L`: Quadrupole length

# Returns
- `mx, my`: Horizontal and vertical transfer matrices
"""
function linear_quad_matrices(K1, L)
  sqrtk = sqrt(abs(K1))
  w = sqrtk*L

  mf = SA[cos(w)        L*sincu(w);
          -sqrtk*sin(w) cos(w)     ]
  
  md = SA[cosh(w)        L*sinhcu(w);
          sqrtk*sinh(w) cosh(w)      ]

  if K1 >= 0
    return mf, md
  else
    return md, mf
  end
end

"""
    linear_thin_quad_matrices(K1L)

Generate transfer matrices for a thin quadrupole.

# Arguments
- `K1L`: Integrated quadrupole strength

# Returns
- `mx, my`: Horizontal and vertical transfer matrices
"""
function linear_thin_quad_matrices(K1L)
  mx = SA[1     0;
          -K1L  1]
  my = SA[1     0;
          K1L   1]

  return mx, my
end 

"""
    linear_solenoid_matrix(Ks, L)

Generate transfer matrix for a solenoid.

# Arguments
- `Ks`: Solenoid strength
- `L`: Solenoid length

# Returns
- 4x4 transfer matrix for coupled horizontal and vertical motion

# Notes
- Based on Bmad manual "Solenoid Tracking" section
"""
function linear_solenoid_matrix(Ks, L)
  s, c = sincos(Ks*L)

  return SA[(1+c)/2     s/Ks       s/2          (1-c)/Ks;
            -Ks*s/4     (1+c)/2    -Ks*(1-c)/4  s/2     ;
            -s/2        -(1-c)/Ks  (1+c)/2      s/Ks    ;
            Ks*(1-c)/4  -s/2       -Ks*s/4      (1+c)/2 ;]
end


"""
    linear_bend_matrices(K0, L, gamma_0, e1, e2)

Generate transfer matrices for a bending magnet.

# Arguments
- `K0`: Bending strength
- `L`: Bend length
- `gamma_0`: Reference gamma
- `e1`: Entrance edge angle (optional)
- `e2`: Exit edge angle (optional)

# Returns
- `mx, my`: Horizontal and vertical transfer matrices
- `r56`: Momentum compaction term
- `d`: Dispersion vector
- `t`: Path length terms
"""
function linear_bend_matrices(K0, L, gamma_0, e1=nothing, e2=nothing)
  theta = K0*L
  s, c = sincos(theta)
  cc = (sincu(theta/2)^2)/2
  sc = sincu(theta)
  mx = SA[c  L*sc; -K0*s  c]
  my = SA[1  L; 0 1]
  r56 = L*(1/gamma_0^2 - theta^2*sincuc(theta))
  d = SA[theta*L*cc, theta*sc, 0, 0]
  t = SA[-theta*sc,  -theta*L*cc, 0, 0]

  if !isnothing(e1) && e1 != 0
    me1 = K0*tan(e1)
    mx = mx*SA[1 0; me1  1]
    my = my*SA[1 0; -me1 1]
    t = SA[t[1]+me1*t[2], t[2], 0, 0]
  end

  if !isnothing(e2) && e2 != 0
    me2 = K0*tan(e2)
    mx = SA[1 0; me2  1]*mx
    my = SA[1 0; -me2 1]*my
    d = SA[d[1], me2*d[1]+d[2], 0, 0]
  end

  return mx, my, r56, d, t
end

end