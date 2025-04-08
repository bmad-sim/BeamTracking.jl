#=

Linear tracking methods expanded around "zero orbit".

=#
# Define the Linear tracking method, and number of rows in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Linear end
MAX_TEMPS(::Linear) = 0

module LinearTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI

const TRACKING_METHOD = Linear

# Drift kernel
@inline function linear_drift!(i, v, v0, work, L, r56)
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = v0[i,XI] + v0[i,PXI] * L
    v[i,YI]  = v0[i,YI] + v0[i,PYI] * L
    v[i,ZI]  = v0[i,ZI] + v0[i,PZI] * r56
    v[i,PXI] = v0[i,PXI]
    v[i,PYI] = v0[i,PYI]
    v[i,PZI] = v0[i,PZI]
  end end
  return v
end


#=
 Generic function for an uncoupled matrix with coasting plane:

[ mx      0       0   d[1:2]]
[ 0       my      0   d[3:4]]
[ t[1:2]  t[3:4]  1   r56   ]

=#
@inline function linear_coast_uncoupled!(i, v, v0, work, mx::AbstractMatrix, my::AbstractMatrix, r56, d::Union{AbstractVector,Nothing}=nothing, t::Union{AbstractVector,Nothing}=nothing)
  @assert size(mx) == (2,2) "Size of matrix mx must be (2,2) for linear_coast_uncoupled!. Received $(size(mx))"
  @assert size(my) == (2,2) "Size of matrix my must be (2,2) for linear_coast_uncoupled!. Received $(size(my))"
  @assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast_uncoupled!. Received $d"
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = mx[1,1] * v0[i,XI] + mx[1,2] * v0[i,PXI] 
    v[i,PXI] = mx[2,1] * v0[i,XI] + mx[2,2] * v0[i,PXI]
    v[i,YI]  = my[1,1] * v0[i,YI] + my[1,2] * v0[i,PYI] 
    v[i,PYI] = my[2,1] * v0[i,YI] + my[2,2] * v0[i,PYI]
    v[i,ZI]  = v0[i,ZI] + r56 * v0[i,PZI]
    v[i,PZI] = v0[i,PZI]
  end end
  if !isnothing(d)
    @inbounds begin @FastGTPSA! begin
      v[i,XI]  += d[XI]  * v0[i,PZI]
      v[i,PXI] += d[PXI] * v0[i,PZI]
      v[i,YI]  += d[YI]  * v0[i,PZI]
      v[i,PYI] += d[PYI] * v0[i,PZI]
    end end
  end
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v0[i,XI] + t[PXI] * v0[i,PXI] + t[YI] * v0[i,YI] + t[PYI] * v0[i,PYI]
    end end
  end
  return v
end

@inline function linear_coast!(i, v, v0, work, mxy::AbstractMatrix, r56, d::Union{AbstractVector,Nothing}=nothing, t::Union{AbstractVector,Nothing}=nothing)
  @assert size(mxy) == (4,4) "Size of matrix mxy must be (4,4) for linear_coast!. Received $(size(mxy))"
  @assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast!. Received $d"
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = mxy[XI, XI] * v0[i,XI] + mxy[XI, PXI] * v0[i,PXI] + mxy[XI, YI] * v0[i,YI] + mxy[XI, PYI] * v0[i,PYI]
    v[i,PXI] = mxy[PXI,XI] * v0[i,XI] + mxy[PXI,PXI] * v0[i,PXI] + mxy[PXI,YI] * v0[i,YI] + mxy[PXI,PYI] * v0[i,PYI]
    v[i,YI]  = mxy[YI, XI] * v0[i,XI] + mxy[YI, PXI] * v0[i,PXI] + mxy[YI, YI] * v0[i,YI] + mxy[YI, PYI] * v0[i,PYI] 
    v[i,PYI] = mxy[PYI,XI] * v0[i,XI] + mxy[PYI,PXI] * v0[i,PXI] + mxy[PYI,YI] * v0[i,YI] + mxy[PYI,PYI] * v0[i,PYI]
    v[i,ZI]  = v0[i,ZI] + r56 * v0[i,PZI]
    v[i,PZI] = v0[i,PZI]
  end end
  if !isnothing(d)
    @inbounds begin @FastGTPSA! begin
      v[i,XI]  += d[XI]  * v0[i,PZI]
      v[i,PXI] += d[PXI] * v0[i,PZI]
      v[i,YI]  += d[YI]  * v0[i,PZI]
      v[i,PYI] += d[PYI] * v0[i,PZI]
    end end
  end
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v0[i,XI] + t[PXI] * v0[i,PXI] + t[YI] * v0[i,YI] + t[PYI] * v0[i,PYI]
    end end
  end
  return v
end

@inline function linear_6D!(i, v, v0, work, m::AbstractMatrix)
  @assert size(m) == (4,4) "Size of matrix m must be (6,6) for linear_6D!. Received $(size(m))"
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = m[XI, XI] * v0[i,XI] + m[XI, PXI] * v0[i,PXI] + m[XI, YI] * v0[i,YI] + m[XI, PYI] * v0[i,PYI] + m[XI, ZI] * v0[i,ZI] + m[XI, PZI] * v0[i,PZI]
    v[i,PXI] = m[PXI,XI] * v0[i,XI] + m[PXI,PXI] * v0[i,PXI] + m[PXI,YI] * v0[i,YI] + m[PXI,PYI] * v0[i,PYI] + m[PXI,ZI] * v0[i,ZI] + m[PXI,PZI] * v0[i,PZI]
    v[i,YI]  = m[YI, XI] * v0[i,XI] + m[YI, PXI] * v0[i,PXI] + m[YI, YI] * v0[i,YI] + m[YI, PYI] * v0[i,PYI] + m[YI, ZI] * v0[i,ZI] + m[YI, PZI] * v0[i,PZI]
    v[i,PYI] = m[PYI,XI] * v0[i,XI] + m[PYI,PXI] * v0[i,PXI] + m[PYI,YI] * v0[i,YI] + m[PYI,PYI] * v0[i,PYI] + m[PYI,ZI] * v0[i,ZI] + m[PYI,PZI] * v0[i,PZI]
    v[i,ZI]  = m[ZI, XI] * v0[i,XI] + m[ZI, PXI] * v0[i,PXI] + m[ZI, YI] * v0[i,YI] + m[ZI, PYI] * v0[i,PYI] + m[ZI, ZI] * v0[i,ZI] + m[ZI, PZI] * v0[i,PZI]
    v[i,PZI] = m[PZI,XI] * v0[i,XI] + m[PZI,PXI] * v0[i,PXI] + m[PZI,YI] * v0[i,YI] + m[PZI,PYI] * v0[i,PYI] + m[PZI,ZI] * v0[i,ZI] + m[PZI,PZI] * v0[i,PZI]
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