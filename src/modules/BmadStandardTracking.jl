struct BmadStandard 
    spin::String
end
MAX_TEMPS(::Magnus) = 8

module BmadStandardTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
const TRACKING_METHOD = BmadStandard


# Pipe/drift - no spin rotation
@makekernel fastgtpsa=true function magnus_drift!(i, b::BunchView, p0c, beta, L)
    # No spin rotation for drift
end

# Solenoid
@makekernel fastgtpsa=true function magnus_solenoid!(i, b::BunchView, ks, γ, G, L)
    v = b.v
    # Constants
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    ι = 2 * γ - 1
    # Check for negligible ks
    if abs(ks) < 1e-20
        return
    end
    pr = sqrt(rel_p^2 - (v[i,PXI] + v[i,YI]*ks/2)^2 - (v[i,PYI] - v[i,XI]*ks/2)^2)
    c = cos(ks*L/pr)
    s = sin(ks*L/pr)
    # Calculate a, b, cc components
    a = 2*pr*(c-1)*v[i,PYI] - 2*pr*s*v[i,PXI] + ks^2*L*v[i,YI]
    a = a - ks*(2*L*v[i,PXI] + (c-1)*pr*v[i,XI] + pr*s*v[i,YI])
    a = a * ξ / (8.0*rel_p^2)
    b = ks*L*(2*v[i,PYI] + ks*v[i,XI])
    b = b + pr*(2*(c-1)*v[i,PXI] + 2*s*v[i,PYI] + c*ks*v[i,YI] - ks*(s*v[i,XI] + v[i,YI]))
    b = -b * ξ / (8.0*rel_p^2)
    cc = 4.0*pr*(v[i,PXI]^2 + v[i,PYI]^2)*s*(1 + G*ι) + ks^3*L*(v[i,XI]^2 + v[i,YI]^2)*(1 + G*ι) -
         ks^2*pr*s*(v[i,XI]^2 + v[i,YI]^2)*(1 + G*ι) + 4.0*ks*((-1+c)*pr*(v[i,PXI]*v[i,XI] + v[i,PYI]*v[i,YI])* 
         (1 + G*ι) + L*(v[i,PXI]^2 + v[i,PYI]^2 + 4.0*rel_p^2 + G*(4.0*rel_p^2 + (v[i,PXI]^2 + v[i,PYI]^2)*ι)))
    cc = cc / (32.0*rel_p^3)
    cc2 = 2*ks^2*L*(v[i,XI]*(-2*v[i,PYI] + 2*v[i,PXI]*s + ks*v[i,XI]) + 2*(v[i,PXI] + v[i,PYI]*s)*v[i,YI] + ks*v[i,YI]^2) + pr*
          (4.0*s*(v[i,PXI]^2 + v[i,PYI]^2) + 4.0*ks*v[i,PYI]*(s*v[i,XI] - 2*v[i,YI]) - 4.0*ks*v[i,PXI]*(2*v[i,XI] + s*v[i,YI]) -
          3.0*ks^2*s*(v[i,XI]^2 + v[i,YI]^2)) + c*ks*(8.0*pr*(v[i,PXI]*v[i,XI] + v[i,PYI]*v[i,YI]) + L*(-4.0*(v[i,PXI]^2 + v[i,PYI]^2) + ks^2*(v[i,XI]^2 + v[i,YI]^2)))
    cc2 = -cc2*pr*ξ^2 / (64.0*rel_p^4)
    cc = cc + cc2
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)
end

# SBend
@makekernel fastgtpsa=true function magnus_sbend!(i, b::BunchView, g, k0, γ0, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    kx = g * k0
    if abs(g) < 1e-20 && abs(k0) < 1e-20
        return
    end
    omegax = sqrt(abs(kx)/rel_p)
    xc = abs(kx) > 0 ? (g*rel_p - k0)/kx : 0.0
    xd = v[i,XI] - xc
    local cx, sx, cx2, sx2, τx
    if kx > 0
        cx = cos(omegax*L)
        sx = sin(omegax*L)/omegax
        cx2 = cos(omegax*L/2)
        sx2 = sin(omegax*L/2)
        τx = -1
    elseif kx < 0
        cx = cosh(omegax*L)
        sx = sinh(omegax*L)/omegax
        cx2 = cosh(omegax*L/2)
        sx2 = sinh(omegax*L/2)
        τx = 1
    else
        cx = 1
        sx = 0.0
        τx = 0.0
    end
    υ = v[i,PYI]^2*ξ - rel_p^2*χ
    pry = sqrt(rel_p^2 - v[i,PYI]^2)
    local a, b, cc
    if g == 0
        a = -0.5*k0*v[i,PYI]*ξ*(L*v[i,PXI] - 0.5*k0*L^2)/(rel_p^2*pry)
        b = -0.5*k0*(L*v[i,PXI]^2 - k0*L^2*v[i,PXI] + k0^2*L^3/3.0)*υ/(2*rel_p^2*pry^3) -
             0.5*k0*L*υ/(rel_p^2*pry)
        cc = -0.5*k0*v[i,PYI]*ξ*L/(rel_p^2)
    elseif k0 == 0
        a = 0.0
        b = -0.5*g*L
        cc = 0.0
    else
        sigma = rel_p*sx*(v[i,XI]-3.0*xc)*xd + L*(-(rel_p*xd*(v[i,XI]-(2+cx)*xc)) + v[i,PXI]*(sx*xc + τx*v[i,PXI]/(rel_p*omegax^2))) -
                v[i,PXI]*(v[i,PXI]*sx + 2*(cx-1)*rel_p*xc)*τx/(rel_p*omegax^2)
        a = -k0*v[i,PYI]*sx2*ξ*(g*v[i,PXI]*sx*omegax + rel_p*(2+g*(v[i,XI]*(1+cx)+xc*(1-cx)))*omegax)*
            (cx2*v[i,PXI] + rel_p*sx2*xd*τx*omegax)/(2*rel_p^3*omegax^2*pry)
        b = -k0*L*υ/(2*rel_p^2*pry)
        b = b + (-g*k0*υ*(sx*xd + L*xc + (cx-1)*v[i,PXI]*τx/(rel_p*omegax^2))/(2*rel_p^2*pry))
        b = b + (-k0*υ*(cx^2*v[i,PXI]*rel_p*xd + v[i,PXI]*rel_p*xd*τx*(-τx + sx^2*omegax^2) + 
                 L*(v[i,PXI]^2 - rel_p^2*xd^2*τx*omegax^2) + cx*sx*(v[i,PXI]^2 + rel_p^2*xd^2*τx*omegax^2))/(8.0*rel_p^2*pry^3))
        b = b + (sigma*g*k0^2*v[i,PYI]^2*ξ^2/(4.0*rel_p^4*pry))
        b = b - 0.5*g*L
        cc = -k0*v[i,PYI]*ξ*(L + g*sx*xd + g*L*xc + (τx*(cx-1)*g*v[i,PXI]/(rel_p*omegax^2)))/(2*rel_p^2)
        cc = cc + (-sigma*g*k0^2*v[i,PYI]*ξ*υ/(4.0*rel_p^4*pry^2))
    end
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)

    # Update coordinates
    BmadStandardexact_sbend!(i, b::BunchView, beta_0, brho_0, hc, b0, e1, e2, Lr)
end

# Quadrupole
@makekernel fastgtpsa=true function magnus_quadrupole!(i, b::BunchView, k1, γ, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    if abs(k1) < 1e-20
        return
    end
    k1_norm = k1 / rel_p
    omega = sqrt(abs(k1_norm))
    s = sin(omega*L)
    sh = sinh(omega*L)
    c = cos(omega*L)
    ch = cosh(omega*L)
    local cx, sx, cy, sy, τy
    if k1_norm > 0
        cx = c
        sx = s/omega
        cy = ch
        sy = sh/omega
        τy = 1
    else
        cy = c
        sy = s/omega
        cx = ch
        sx = sh/omega
        τy = -1
    end
    a = 0.5*χ*k1_norm*(sy*v[i,YI] + τy*(cy-1)/(omega^2)*v[i,PYI]/rel_p)
    b = 0.5*χ*k1_norm*(sx*v[i,XI] + τy*(1-cx)/(omega^2)*v[i,PXI]/rel_p)
    cc = v[i,PXI]*sx*(v[i,PYI]*sy + cy*rel_p*v[i,YI]) + rel_p*v[i,XI]*(cx*v[i,PYI]*sy + rel_p*(cx*cy-1)*v[i,YI])
    cc = cc * (-0.5*ξ*k1_norm/rel_p^2)
    cc = cc + 0.25*χ^2/rel_p^2*((cx-cy)*v[i,PXI]*v[i,PYI] + τy*omega^2*
         (v[i,PXI]*v[i,PYI]*sx*sy + v[i,PXI]*rel_p*(cy*sx-sy)*v[i,YI] + rel_p*v[i,XI]*(v[i,PYI]*(cx*sy-sx) + (cx*cy-1)*rel_p*v[i,YI])))
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)
end

# Sextupole
@makekernel fastgtpsa=true function magnus_sextupole!(i, b::BunchView, k2, γ, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)

    # First kick of a kick-drift-kick split
    px1 = v[i,PXI] + 0.25*k2*L*(v[i,YI]^2 - v[i,XI]^2)
    py1 = v[i,PYI] + 0.5*k2*L*v[i,XI]*v[i,YI]
    pl = sqrt(rel_p^2 - px1^2 - py1^2)
    a = px1*ξ*(L^2*(-3.0*px1^2*py1 + py1^3) - 
        3.0*L*pl*(2*px1*py1*v[i,XI] + v[i,YI]*(px1^2 - py1^2)) -
        3.0*pl^2*(2*px1*v[i,XI]*v[i,YI] + py1*(v[i,XI]^2 - v[i,YI]^2)))
    a = a + rel_p^2*χ*(L*py1*(2*L*px1 + 3.0*pl*v[i,XI]) + 3.0*pl*(L*px1 + 2*pl*v[i,XI])*v[i,YI])
    a = a * k2*L/(12.0*pl^3*rel_p^2)
    b = py1*ξ*(L^2*(-3.0*px1^2*py1 + py1^3) -
        3.0*L*pl*(2*px1*py1*v[i,XI] + v[i,YI]*(px1^2 - py1^2)) -
        3.0*pl^2*(2*px1*v[i,XI]*v[i,YI] + py1*(v[i,XI]^2 - v[i,YI]^2)))
    b = b + rel_p^2*χ*(L^2*(px1^2 - py1^2) + 3.0*pl^2*(v[i,XI]^2 - v[i,YI]^2) + 3.0*L*pl*(px1*v[i,XI] - py1*v[i,YI]))
    b = b * k2*L/(12.0*pl^3*rel_p^2)
    cc = L^2*(-3.0*px1^2*py1 + py1^3) - 3.0*L*pl*(2*px1*py1*v[i,XI] + v[i,YI]*(px1^2 - py1^2)) -
         3.0*pl^2*(2*px1*v[i,XI]*v[i,YI] + py1*(v[i,XI]^2 - v[i,YI]^2))
    cc = cc * k2*L*ξ/(12.0*pl^2*rel_p^2)
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)
end

end