#=

Linear tracking methods expanded around "zero orbit".

=#
# Define the Linear tracking method, and number of rows in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Linear end
MAX_TEMPS(::Linear) = 5

module LinearTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
const TRACKING_METHOD = Linear

# Maybe get rid of inline here and put in function-wise launch! ?
# Drift kernel
@makekernel function linear_drift!(i, v, L, r56)
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
@makekernel function linear_coast_uncoupled!(i, v, mx::AbstractMatrix, my::AbstractMatrix, r56, d::Union{AbstractArray,Nothing}, t::Union{AbstractArray,Nothing})
  #@assert size(work, 2) >= 1 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 1) for linear_coast_uncoupled!"
  #@assert size(mx) == (2,2) "Size of matrix mx must be (2,2) for linear_coast_uncoupled!. Received $(size(mx))"
  #@assert size(my) == (2,2) "Size of matrix my must be (2,2) for linear_coast_uncoupled!. Received $(size(my))"
  #@assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast_uncoupled!. Received $d"
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
    end end
  end
  w1 = v[i,XI]
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = mx[1,1] * v[i,XI]   + mx[1,2] * v[i,PXI] 
    v[i,PXI] = mx[2,1] * w1 + mx[2,2] * v[i,PXI]
    w1 = v[i,YI]
    v[i,YI]  = my[1,1] * v[i,YI]   + my[1,2] * v[i,PYI] 
    v[i,PYI] = my[2,1] * w1 + my[2,2] * v[i,PYI]
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

@makekernel function linear_coast!(i, v, mxy::AbstractMatrix, r56, d::Union{AbstractArray,Nothing}, t::Union{AbstractArray,Nothing})
  #@assert size(work, 2) >= 3 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 3) for linear_coast!"
  #@assert size(mxy) == (4,4) "Size of matrix mxy must be (4,4) for linear_coast!. Received $(size(mxy))"
  #@assert isnothing(d) || length(d) == 4 "The dispersion vector d must be either `nothing` or of length 4 for linear_coast!. Received $d"
  if !isnothing(t)
    @inbounds begin @FastGTPSA! begin
      v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
    end end
  end
  w1= v[i,XI]
  w2= v[i,PXI]
  w3= v[i,YI]
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = mxy[XI, XI] * v[i,XI]   + mxy[XI, PXI] * v[i,PXI]  + mxy[XI, YI] * v[i,YI]   + mxy[XI, PYI] * v[i,PYI]
    v[i,PXI] = mxy[PXI,XI] * w1 + mxy[PXI,PXI] * v[i,PXI]  + mxy[PXI,YI] * v[i,YI]   + mxy[PXI,PYI] * v[i,PYI]
    v[i,YI]  = mxy[YI, XI] * w1 + mxy[YI, PXI] * w2 + mxy[YI, YI] * v[i,YI]   + mxy[YI, PYI] * v[i,PYI] 
    v[i,PYI] = mxy[PYI,XI] * w1 + mxy[PYI,PXI] * w2 + mxy[PYI,YI] * w3 + mxy[PYI,PYI] * v[i,PYI]
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

@makekernel function linear_6D!(i, v, m::AbstractMatrix)
  #@assert size(work, 2) >= 5 && size(work, 1) >= size(v, 1) "Size of work matrix must be at least ($(size(v, 1)), 5) for linear_6D!"
  #@assert size(m) == (6,6) "Size of matrix m must be (6,6) for linear_6D!. Received $(size(m))"
  w1= v[i,XI]
  w2= v[i,PXI]
  w3= v[i,YI]
  w4= v[i,PYI]
  w5= v[i,ZI]
  @inbounds begin @FastGTPSA! begin
    v[i,XI]  = m[XI, XI] * v[i,XI]   + m[XI, PXI] * v[i,PXI]  + m[XI, YI] * v[i,YI]   + m[XI, PYI] * v[i,PYI]  + m[XI, ZI] * v[i,ZI]   + m[XI, PZI] * v[i,PZI]
    v[i,PXI] = m[PXI,XI] * w1 + m[PXI,PXI] * v[i,PXI]  + m[PXI,YI] * v[i,YI]   + m[PXI,PYI] * v[i,PYI]  + m[PXI,ZI] * v[i,ZI]   + m[PXI,PZI] * v[i,PZI]
    v[i,YI]  = m[YI, XI] * w1 + m[YI, PXI] * w2 + m[YI, YI] * v[i,YI]   + m[YI, PYI] * v[i,PYI]  + m[YI, ZI] * v[i,ZI]   + m[YI, PZI] * v[i,PZI]
    v[i,PYI] = m[PYI,XI] * w1 + m[PYI,PXI] * w2 + m[PYI,YI] * w3 + m[PYI,PYI] * v[i,PYI]  + m[PYI,ZI] * v[i,ZI]   + m[PYI,PZI] * v[i,PZI]
    v[i,ZI]  = m[ZI, XI] * w1 + m[ZI, PXI] * w2 + m[ZI, YI] * w3 + m[ZI, PYI] * w4 + m[ZI, ZI] * v[i,ZI]   + m[ZI, PZI] * v[i,PZI]
    v[i,PZI] = m[PZI,XI] * w1 + m[PZI,PXI] * w2 + m[PZI,YI] * w3 + m[PZI,PYI] * w4 + m[PZI,ZI] * w5 + m[PZI,PZI] * v[i,PZI]
  end end
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

function linear_thin_quad_matrices(K1L)
  mx = SA[1     0;
          -K1L  1]
  my = SA[1     0;
          K1L   1]

  return mx, my
end 

# From the Bmad manual "Solenoid Tracking" section, linearized
function linear_solenoid_matrix(Ks, L)
  s, c = sincos(Ks*L)

  return SA[(1+c)/2     s/Ks       s/2          (1-c)/Ks;
            -Ks*s/4     (1+c)/2    -Ks*(1-c)/4  s/2     ;
            -s/2        -(1-c)/Ks  (1+c)/2      s/Ks    ;
            Ks*(1-c)/4  -s/2       -Ks*s/4      (1+c)/2 ;]
end


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
