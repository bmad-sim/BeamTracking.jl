#=

Linear tracking methods expanded around "zero orbit".

=#
# Define the Linear tracking method, and number of rows in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Linear end

MAX_TEMPS(::Linear) = 1

module LinearTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI

const TRACKING_METHOD = Linear

# Maybe get rid of inline here and put in function-wise launch! ?
# Drift kernel
@inline function linear_drift!(i, v, work, L, r56)
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
@inline function linear_coast_uncoupled!(i, v, work, mx::AbstractMatrix, my::AbstractMatrix, r56, d::Union{AbstractVector,Nothing}=nothing, t::Union{AbstractVector,Nothing}=nothing)
  @assert size(work, 2) >= 1 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 1) for linear_coast_uncoupled!"
  @assert size(mx) == (2,2) "Size of matrix mx must be (2,2) for linear_coast_uncoupled!. Received $(size(mx))"
  @assert size(my) == (2,2) "Size of matrix my must be (2,2) for linear_coast_uncoupled!. Received $(size(my))"
  @assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast_uncoupled!. Received $d"
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

@inline function linear_coast!(i, v, work, mxy::AbstractMatrix, r56, d::Union{AbstractVector,Nothing}=nothing, t::Union{AbstractVector,Nothing}=nothing)
  @assert size(work, 2) >= 3 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 3) for linear_coast!"
  @assert size(mxy) == (4,4) "Size of matrix mxy must be (4,4) for linear_coast!. Received $(size(mxy))"
  @assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast!. Received $d"
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

@inline function linear_6D!(i, v, work, m::AbstractMatrix)
  @assert size(work, 2) >= 5 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 5) for linear_6D!"
  @assert size(m) == (6,6) "Size of matrix m must be (6,6) for linear_6D!. Received $(size(m))"
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
  return v
end

# Utility functions to create a linear matrix
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
#=
# Quadrupole kernel
@inline function linear_quad!(i, v, work, fqi, fpi, dqi, dpi, sqrtk, L, gamma_0)
  @assert all(t->t<=4 && t>=1, (fqi, fpi, dqi, dpi)) "Invalid focus/defocus indices for quadrupole"
  @assert size(work, 2) >= 1 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 1) for linear_quad!"
  @inbounds begin
    @FastGTPSA! work[i,1] = 0 + v[i,fqi]
    @FastGTPSA! v[i,fqi]  = cos(sqrtk*L)*v[i,fqi] + L*sincu(sqrtk*L)*v[i,fpi]
    @FastGTPSA! v[i,fpi]  = -sqrtk*sin(sqrtk*L)*work[i,1] + cos(sqrtk*L)*v[i,fpi]
    @FastGTPSA! work[i,1] = 0 + v[i,dqi]
    @FastGTPSA! v[i,dqi]  = cosh(sqrtk*L)*v[i,dqi] + L*sinhcu(sqrtk*L)*v[i,dpi]
    @FastGTPSA! v[i,dpi]  = sqrtk*sinh(sqrtk*L)*work[i,1] + cosh(sqrtk*L)*v[i,dpi]
    @FastGTPSA! v[i,ZI]   = v[i,ZI] + v[i,PZI] * L/gamma_0^2
  end 
  return v
end
=#


end