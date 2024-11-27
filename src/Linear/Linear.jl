module Linear
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

"""
    struct Linear.Drift{T}
  
## Fields
• `L` -- Drift length / m \\
"""
Base.@kwdef struct Drift{T}
  L::T  # drift length / m
end

"""
    struct Linear.Quadrupole{T}

## Fields
• `L`  -- Quadrupole length / m \\
• `Bn1` -- Quadrupole gradient / (T·m^-1)
"""
Base.@kwdef struct Quadrupole{T}
  L::T   # quadrupole length / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end


Base.@kwdef struct SBend{T}
  L::T   # arc length / m
  B0::T  # field strength / T
  g::T   # coordinate system curvature through element / m^-1
  e1::T  # edge 1 angle / rad
  e2::T  # edge 2 angle / rad
end

Base.@kwdef struct Combined{T}
  L::T 
  B0::T 
  Bn1::T
  g::T
  e1::T
  e2::T
end

Base.@kwdef struct Solenoid{T}
  L::T
  Bs::T
end


function track!(bunch::Bunch, ele::Linear.Drift; work=get_work(bunch, Val{0}()))
  L = ele.L
  v = bunch.v
  
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)

  @FastGTPSA! begin
  @. v.x  = v.x + L*v.px
  @. v.y  = v.y + L*v.py
  @. v.z  = v.z + L/gamma_ref^2*v.pz
  end

  # Spin unchanged
  
  return bunch
end


function track!(bunch::Bunch, ele::Linear.Quadrupole; work=get_work(bunch, Val{1}()))
  v = bunch.v
  L = ele.L

  Kn1 = ele.Bn1 / brho(massof(bunch.species), bunch.beta_gamma_ref, chargeof(bunch.species))
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)

  if Kn1 >= 0
    k = sqrt(Kn1)
    kL = k*L
    sx  = sin(kL)
    cx  = cos(kL)
    sxc = sincu(kL)
    sy  = sinh(kL)
    cy  = cosh(kL)
    syc = sinhcu(kL)
    sgn = 1
  else
    k = sqrt(-Kn1)
    kL = k*L
    sx  = sinh(kL)
    cx  = cosh(kL)
    sxc = sinhcu(kL)
    sy  = sin(kL)
    cy  = cos(kL)
    syc = sincu(kL)
    sgn = -1
  end

  
  # Note: 0+work[1] is workaround for GTPSA @FastGTPSA! bug
  @FastGTPSA! begin
  @. work[1]  = cx*v.x + sxc*L*v.px
  @. v.px     = -sgn*k*sx*v.x + cx*v.px
  @. v.x      = 0+work[1] 

  @. work[1]  = cy*v.y + syc*L*v.py
  @. v.py     = sgn*k*sy*v.y + cy*v.py
  @. v.y      = 0+work[1] 

  @. v.z      = v.z + L/gamma_ref^2 * v.pz
  end 

  return bunch
end 

"""
Routine to linearly tracking through a Sector Bending Magnet
"""
function track!(bunch::Bunch, ele::Linear.SBend; work=getwork(bunch,Val{1}()))
  v = bunch.v
  L = ele.L
  e1 = ele.e1
  e2 = ele.e2
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)


  #curvature of B field -- confusion: g =  1/rho = qB0/p? 
  gtot = ele.B0 / brho(massof(bunch.species),bunch.beta_gamma_ref,chargeof(bunch.species))

  pd = g * gtot
  K = sqrt(pd)
  kL = K * L
  diff =  g - gtot
  
  co = cos(kL)
  si = sin(kL)
  tan1 = tan(e1)
  tan2 = tan(e2)
  
  @FastGTPSA! begin
    @. m21 = K * ((tan1 + tan2) * co + (tan1 * tan2 - 1) * si)
    @. m26 = (g + diff)* si / K + K/gtot *(1 - co)* tan2

    @. m51 = - g/K * si + (diff)*l^2/2 * (1-co*si)/kL - (K*(1-co)/gtot - (diff)*si^2/K/2) * tan1
    @. m52 = (co - 1) / gtot - (diff) * (si^2 / pd / 2)
    @. m56 = - g / gtot * L + L/gamma_ref^2 + g^2 / K^3 * si - (diff)^2 * L * (1 - co * si / kL) / gtot / 2
    @. cons =  - diff / gtot * L + g * (diff) * si / K^3 - (diff)^2 * L * (1 - co * si / kL) /pd / 4

    @. work[1] = (1 - kL * tan1) * v.y + L * v.py #new y
    @. v.py = - K * (tan2 + tan1 * (1 - kL *tan2)) * v.y + (1 - kL * tan2) * v.py 
    @. v.y = 0 + work[1]

    @. work[1] = (co + si * tan1) * v.x + si / K * v.px + (1 - co)/gtot * v.pz + (diff)*(1 - co)/pd #new x
    @. v.px = m21 * v.x + (co + si * tan2) *  v.px + m26 * v.pz + (diff) * si / K
    @. v.z = m51 * v.x + m52 * v.px + v.z + m56 * v.pz + cons

    @. v.x = 0 + work[1]
 
  # @. work[1] = (1 - gL * tan1) * v.y + L * v.py  #new y 
  # @. v.py = - g * (tan2 + tan1 * (1 - gL * tan2)) * v.y + (1 - gL * tan2) * v.py
  # @. v.y = 0 + work[1]

  # @. work[1] = (co + si * tan1) * v.x + si / g * v.px + (1 - co) / g * v.pz # new x
  # @. v.px = g * ((tan1 + tan2) * co + (tan1 * tan2 - 1) * si) * v.x + (co + si * tan2) * v.px + (si + (1 - co) * tan2) * pz
  # @. v.z = - (si + (co - 1) * tan1) * v.x + (co - 1) / g * v.px + v.z + (si - gL) / g * v.pz

  # @.v.x = 0 + work[1]

  end 

  return bunch
end 


"""
Routine to linearly tracking through a Combined Magnet
"""
# #function track!(bunch::Bunch, ele::Linear.Combined; work = )
#   v = bunch.v
#   L = ele.L
  
#   brho = brho(massof(bunch.species),bunch.beta_gamma_ref,chargeof(bunch.species))
#   ka = ele.B0 / brho #curvature
#   k1 =  ele.Bn1 / brho #quad strength

#   K = ka^2 + k1

#   sqk1 = sqrt(abs(k1)) 
#   sqK = sqrt(abs(K))

#   Kl = sqK * L
#   kl  = sqk1 * L

#   @FastGTPSA! begin
#   @. zf.x  = zi.x * cos(Kl) +  zi.px * sin(Kl) / sqK + zi.pz * (1 - cos(Kl)) * ka / K
#   @. zf.px = - zi.x * sqK * sin(Kl) + zi.px * cos(Kl) + zi.pz * sin(Kl) * ka / sqK
#   @. zf.y =  zi.y * cosh(kl) + zi.py * sinh(kl) / sqk
#   @. zf.py = zi.y * sqK * sinh(kl) + zi.py * cosh(kl)
#   @. zf.z = - zi.x * sin(Kl) * ka / sqK + zi.px * (cos(Kl)-1) * ka / K + zi.z + zi.pz * (sin(Kl) / sqK - l) *  ka^2 / K
#   @. zf.pz = zi.pz
#   end 
# #end 
# #end



function track!(bunch::Bunch, ele::Linear.Solenoid; work=getwork(bunch,Val{3}()))
  v = bunch.v
  L = ele.L
  S = ele.Bs / brho(massof(bunch.species),bunch.beta_gamma_ref,chargeof(bunch.species))
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)
  phi = S * L / 2

  co = cos(phi)
  si = sin(phi)
  si2 = sin(2*phi)

  @FastGTPSA!  begin
  @. work[1] = co^2 * v.x + si2 / S * v.px + si2 / 2 * v.y + si^2 * 2 / S * v.py #new x
  @. work[2] = - si2 * S / 4 * v.x + co^2 * v.px - si^2 * S / 2 * v.y + si2 / 2 * v.py #new px 
  @. work[3] = - si2 / 2 * v.x - si^2 * 2 / S * v.px + co^2 * v.y + si2 / S * v.py
  @. v.py = si^2 * S / 2 * v.x - si2 / 2 * v.px - si2 * S / 4 * v.y + co^2 * v.py 

  @. v.x = 0 + work[1]
  @. v.px = 0 + work[2]
  @. v.y = 0 + work[3]
  @. v.z  = v.z + L/gamma_ref^2 * v.pz
  end

  return bunch
end

end