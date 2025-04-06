module Linear
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

"""
    struct Linear.Drift{T}
  
## Fields
• `L` -- Drift L / m \\
"""
Base.@kwdef struct Drift{T}
  L::T  # drift L / m
end

"""
    struct Linear.Quadrupole{T}

## Fields
• `L`  -- Quadrupole L / m \\
• `Bn1` -- Quadrupole gradient / (T·m^-1)
"""
Base.@kwdef struct Quadrupole{T}
  L::T   # quadrupole L / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end


Base.@kwdef struct SBend{T}
  L::T   # arc L / m
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
  cons = - dg / gtot * L + g * dg * si / K^3 - dg^2 * L * (1 - co * si / kL) /pd / 4
  
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
  

  if k1 >= 0
    wy = sqrt(k1)
    wyL = wy * L 
    cy  = cosh(wyL)
    syc = sinhcu(wyL) * L
    sgny = 1
  else
    wy = sqrt(-k1)
    wyL = wy * L 
    cy  = cos(wyL)
    syc = sincu(wyL) * L 
    sgny = - 1
  end

  if kx >= 0 
    wx = sqrt(kx)
    wxL = wx * L 
    cx  = cos(wxL)
    sxc = sincu(wxL) * L 
    sgnx = - 1
  else
    wx = sqrt(-kx)
    wxL = wx * L 
    cx  = cosh(wxL)
    sxc = sinhcu(wxL) * L
    sgnx = 1
  end

  if abs(kx)<1e-10
    z2 = g * L * L / 2
  else 
    z2 = sgnx * g * (1 - cx) / (wx * wx)
  end
   
  
  x_c = (g - gtot)/kx
  dx_c = g/kx

  dom_x = - wx /2
  dom_xx = -1/2
  
  dc_x = sgnx * sxc * wx * dom_x * L
  ds_x = (cx * L - sxc) * dom_xx
  dcs_x = cx * ds_x + dc_x * sxc

  z0  = -g * x_c * L
  z1  = -g * sxc
  z11 = sgnx * wx * wx * (L - cx * sxc) / 4
  z12 = -sgnx * wx * wx * sxc * sxc / 2


  ht_x1 = gtot * tan(e1)
  ht_y1 = -ht_x1

  ht_x2 = gtot * tan(e2)
  ht_y2 = - ht_x2

  @FastGTPSA! begin
  
    @. work[1] = cy * v.y + syc * v.py
    @. v.py = sgny * wy * wy * syc * v.y + cy * v.py + (ht_y1 + ht_y2) * work[1]
    @. v.y = 0 + work[1]
    @. v.z = (z1 - x_c * 2 * z11) * v.x + (z2 - x_c * z12) * v.px + v.z + (- dx_c * (z1 - x_c * 2 * z11) - g * L * dx_c + x_c * g * ds_x - (z11 + sgnx * wx * wx * dcs_x / 4) * x_c * x_c + L/gamma_ref^2) * v.pz
    @. work[1] = cx * v.x + sxc * v.px + (dx_c * (1 - cx) - dc_x * x_c) * v.pz #new x
    @. v.px = sgnx * wx * wx * sxc * v.x + cx * v.px - sgnx * wx * (2 * dom_x * sxc * x_c + wx * sxc * x_c + wx * ds_x * x_c + dx_c * wx * sxc) * v.pz + (ht_x1 + ht_x2) * work[1]
    @. v.x = work[1]
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
