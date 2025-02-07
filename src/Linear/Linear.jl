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
  g = ele.g

  gamma_ref = sr_gamma(bunch.beta_gamma_ref)


  #curvature of B field --  gtot =  1/rho, actual field 
  gtot = ele.B0 / brho(massof(bunch.species),bunch.beta_gamma_ref,chargeof(bunch.species))

  pd = g * gtot
  K = sqrt(pd)
  kL = K * L
  dg = g - gtot
  
  co = cos(kL)
  si = sin(kL)
  tan1 = tan(e1)
  tan2 = tan(e2)

  m21 = K * ((tan1 + tan2) * co + (tan1 * tan2 - 1) * si)
  m26 = (g + dg)* si / K + K/gtot *(1 - co)* tan2
  m51 = - g/K * si + dg*L^2/2 * (1-co*si)/kL - (K*(1-co)/gtot - dg*si^2/K/2) * tan1
  m52 = (co - 1) / gtot - dg * si^2 / pd / 2
  m56 = - g / gtot * L + L/gamma_ref^2 + g^2 / K^3 * si - dg^2 * L * (1 - co * si / kL) / gtot / 2
  cons =  - dg / gtot * L + g * dg * si / K^3 - dg^2 * L * (1 - co * si / kL) /pd / 4
  
  @FastGTPSA! begin
    @. work[1] = (1 - kL * tan1) * v.y + L * v.py #new y
    @. v.py = - K * (tan2 + tan1 * (1 - kL *tan2)) * v.y + (1 - kL * tan2) * v.py 
    @. v.y = 0 + work[1]

    @. work[1] = (co + si * tan1) * v.x + si / K * v.px + (1 - co)/gtot * v.pz + (dg)*(1 - co)/pd #new x
    @. v.z = m51 * v.x + m52 * v.px + v.z + m56 * v.pz + cons
    @. v.px = m21 * v.x + (co + si * tan2) * v.px + m26 * v.pz + dg * si / K
    @. v.x = 0 + work[1]


  end 
  return bunch
end 


"""
Routine to linearly tracking through a Combined Magnet
"""
function track!(bunch::Bunch, ele::Linear.Combined; work=get_work(bunch, Val{1}()))
  v = bunch.v
  L = ele.L
  e1 = ele.e1
  e2 = ele.e2
  g = ele.g

  gamma_ref = sr_gamma(bunch.beta_gamma_ref)
  br = brho(massof(bunch.species),bunch.beta_gamma_ref,chargeof(bunch.species))
  
  gtot = ele.B0 / br #curvature of B field
  k1 =  ele.Bn1 / br #quad strength 
  kx = k1 + g * gtot
  wx2 = abs(kx)
  tan1 = tan(e1)
  tan2 = tan(e2)
 

  if k1 >= 0
    wy = sqrt(k1)
    wyL = wy * L 
    cy  = cosh(wyL)
    sy = sinh(wyL) / wy #if not needed, delete after all sy 
    syc = sinhcu(wyL) * L
    sgny = 1
  else
    wy = sqrt(-k1)
    wyL = wy * L 
    cy  = cos(wyL)
    sy = sin(wyL) / wy
    syc = sincu(wyL) * L 
    sgny = - 1
  end

  if kx >= 0 
    wx = sqrt(kx)
    wxL = wx * L 
    cx  = cos(wxL)
    sx = sin(wxL) / wx
    sxc = sincu(wxL) * L 
    sgnx = - 1
  else
    wx = sqrt(-kx)
    wxL = wx * L 
    cx  = cosh(wxL)
    sxc = sinhcu(wxL) * L
    sgnx = 1
  end
   
  m51 = - g * sxc - (g-gtot) * (L - cx * sxc) * sgnx * wx2 / 2 / kx + gtot * (g * (1-cx) * sgnx / wx2 + (g-gtot) * sxc * sxc * sgnx * wx2 / 2 / kx) * tan1 
  m52 = g*(1-cx)*sgnx/wx/wx + (g-gtot)* sxc * sxc * sgnx * wx2 / 2/ kx
  m56 = g*g*(-L + sxc)/kx + g * (g-gtot)*(L - cx * sxc)* sgnx * wx2 /2/kx/kx + L/gamma_ref^2 
  cons = (g-gtot) /kx * ( - L  + sxc * g + sgnx * wx2 * (L - cx * sxc) * (g - gtot) / kx /4) 

  @FastGTPSA! begin
    @. work[1] = (cy - gtot * syc * tan1) * v.y + syc * v.py #new y 
    @. v.py = (syc * sgny * wy * wy - gtot * cy * tan2 - gtot * tan1 * (cy - gtot * syc * tan2)) * v.y + (cy - gtot * syc * tan2) * v.py 
    @. v.y = 0 + work[1]

    @. v.z = m51 * v.x + m52 * v.px + v.z + m56 * v.pz + cons
    @. work[1] = (cx + g * sxc * tan1) * v.x + sxc * v.px + g * (1 - cx)/kx * v.pz + (1 - cx) * (g-gtot)/kx #new x 
    @. v.px = (sgnx * sxc * wx2 + gtot * cx * tan2 + gtot * tan1 * (cx + gtot * sx * tan2)) * v.x + (cx + gtot * sxc * tan2) * v.px + (-sgnx * wx2 * sxc * (g/kx + (g - gtot)/kx) + g * gtot * (1 - cx) * tan2 / kx) * v. pz - sgnx * wx2 * sxc * (g - gtot) / kx
    @. v.x = 0 + work[1] 
  
  end
  return Bunch
  

end




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