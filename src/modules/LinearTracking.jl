#=

Linear tracking methods expanded around "zero orbit".

=#
# Define the Linear tracking method, and number of rows in the work matrix 
# (equal to number of temporaries needed for a single particle)
struct Linear end

module LinearTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
const TRACKING_METHOD = Linear

# Maybe get rid of inline here and put in function-wise launch! ?
# Drift kernel
@makekernel fastgtpsa=true function linear_drift!(i, b::BunchView, L, r56)
  v = b.v
  v[i,XI] += v[i,PXI] * L
  v[i,YI] += v[i,PYI] * L
  v[i,ZI] += v[i,PZI] * r56
end

#=
 Generic function for an uncoupled matrix with coasting plane:

[ mx      0       0   d[1:2]]
[ 0       my      0   d[3:4]]
[ t[1:2]  t[3:4]  1   r56   ]


=#

@makekernel fastgtpsa=true function linear_coast_uncoupled!(i, b::BunchView, mx::StaticMatrix{2,2}, my::StaticMatrix{2,2}, r56, d::Union{StaticVector{4},Nothing}, t::Union{StaticVector{4},Nothing})
  v = b.v
  if !isnothing(t)
    v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
  end
  old_q = v[i,XI]
  v[i,XI]  = mx[1,1] * v[i,XI] + mx[1,2] * v[i,PXI] 
  v[i,PXI] = mx[2,1] * old_q  + mx[2,2] * v[i,PXI]
  old_q = v[i,YI]
  v[i,YI]  = my[1,1] * v[i,YI] + my[1,2] * v[i,PYI] 
  v[i,PYI] = my[2,1] * old_q  + my[2,2] * v[i,PYI]
  v[i,ZI] += r56 * v[i,PZI]
  if !isnothing(d)
    v[i,XI]  += d[XI]  * v[i,PZI]
    v[i,PXI] += d[PXI] * v[i,PZI]
    v[i,YI]  += d[YI]  * v[i,PZI]
    v[i,PYI] += d[PYI] * v[i,PZI]
  end
end

@makekernel fastgtpsa=true function linear_coast!(i, b::BunchView, mxy::StaticMatrix{4,4}, r56, d::Union{StaticVector{4},Nothing}, t::Union{StaticVector{4},Nothing})
  v = b.v
  if !isnothing(t)
    v[i,ZI] += t[XI] * v[i,XI] + t[PXI] * v[i,PXI] + t[YI] * v[i,YI] + t[PYI] * v[i,PYI]
  end
  old_x  = v[i,XI]
  old_px = v[i,PXI]
  old_y  = v[i,YI]
  v[i,XI]  = mxy[XI, XI] * v[i,XI] + mxy[XI, PXI] * v[i,PXI] + mxy[XI, YI] * v[i,YI] + mxy[XI, PYI] * v[i,PYI]
  v[i,PXI] = mxy[PXI,XI] * old_x  + mxy[PXI,PXI] * v[i,PXI] + mxy[PXI,YI] * v[i,YI] + mxy[PXI,PYI] * v[i,PYI]
  v[i,YI]  = mxy[YI, XI] * old_x  + mxy[YI, PXI] * old_px  + mxy[YI, YI] * v[i,YI] + mxy[YI, PYI] * v[i,PYI] 
  v[i,PYI] = mxy[PYI,XI] * old_x  + mxy[PYI,PXI] * old_px  + mxy[PYI,YI] * old_y  + mxy[PYI,PYI] * v[i,PYI]
  v[i,ZI] += r56 * v[i,PZI]
  if !isnothing(d)
    v[i,XI]  += d[XI]  * v[i,PZI]
    v[i,PXI] += d[PXI] * v[i,PZI]
    v[i,YI]  += d[YI]  * v[i,PZI]
    v[i,PYI] += d[PYI] * v[i,PZI]
  end
end

@makekernel fastgtpsa=true function linear_6D!(i, b::BunchView, m::StaticMatrix{6,6})
  v = b.v
  old_x  = v[i,XI]
  old_px = v[i,PXI]
  old_y  = v[i,YI]
  old_py = v[i,PYI]
  old_z  = v[i,ZI]
  v[i,XI]  = m[XI, XI] * v[i,XI] + m[XI, PXI] * v[i,PXI] + m[XI, YI] * v[i,YI] + m[XI, PYI] * v[i,PYI] + m[XI, ZI] * v[i,ZI] + m[XI, PZI] * v[i,PZI]
  v[i,PXI] = m[PXI,XI] * old_x  + m[PXI,PXI] * v[i,PXI] + m[PXI,YI] * v[i,YI] + m[PXI,PYI] * v[i,PYI] + m[PXI,ZI] * v[i,ZI] + m[PXI,PZI] * v[i,PZI]
  v[i,YI]  = m[YI, XI] * old_x  + m[YI, PXI] * old_px  + m[YI, YI] * v[i,YI] + m[YI, PYI] * v[i,PYI] + m[YI, ZI] * v[i,ZI] + m[YI, PZI] * v[i,PZI]
  v[i,PYI] = m[PYI,XI] * old_x  + m[PYI,PXI] * old_px  + m[PYI,YI] * old_y  + m[PYI,PYI] * v[i,PYI] + m[PYI,ZI] * v[i,ZI] + m[PYI,PZI] * v[i,PZI]
  v[i,ZI]  = m[ZI, XI] * old_x  + m[ZI, PXI] * old_px  + m[ZI, YI] * old_y  + m[ZI, PYI] * old_py  + m[ZI, ZI] * v[i,ZI] + m[ZI, PZI] * v[i,PZI]
  v[i,PZI] = m[PZI,XI] * old_x  + m[PZI,PXI] * old_px  + m[PZI,YI] * old_y  + m[PZI,PYI] * old_py  + m[PZI,ZI] * old_z  + m[PZI,PZI] * v[i,PZI]
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