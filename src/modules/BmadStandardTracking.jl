struct BmadStandard 
    spin::String
end
MAX_TEMPS(::BmadStandard) = 8

module BmadStandardTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
const TRACKING_METHOD = BmadStandard


# Drift - no spin rotation
@makekernel fastgtpsa=true function magnus_drift!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
    ExactTracking.exact_drift!(i, b::BunchView, beta_0, gamsqr_0, tilde_m, L)
end

# Solenoid -- second order Magnus expansion
@makekernel fastgtpsa=true function magnus_solenoid!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, G, L)
    v = b.v
    # Constants
    rel_p = 1 + v[i,PZI]
    γ = sqrt(gamsqr_0)
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

    # Update coordinates
    ExactTracking.exact_solenoid!(i, b::BunchView, ks, beta_0, gamsqr_0, tilde_m, L)
end

# SBend
@makekernel fastgtpsa=true function magnus_sbend!(i, b::BunchView, g, k0, γ0, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    kx = g * k0

    ωx = sqrt(abs(kx)/rel_p)
    xc = abs(kx) > 0 ? (g*rel_p - k0)/kx : 0.0
    xd = v[i,XI] - xc
    if kx > 0
        cx = cos(ωx*L)
        sx = sin(ωx*L)/ωx
        cx2 = cos(ωx*L/2)
        sx2 = sin(ωx*L/2)
        τx = -1
    elseif kx < 0
        cx = cosh(ωx*L)
        sx = sinh(ωx*L)/ωx
        cx2 = cosh(ωx*L/2)
        sx2 = sinh(ωx*L/2)
        τx = 1
    else
        cx = 1
        sx = 0.0
        τx = 0.0
    end
    υ = v[i,PYI]^2*ξ - rel_p^2*χ
    pry = sqrt(rel_p^2 - v[i,PYI]^2)
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
        sigma = rel_p*sx*(v[i,XI]-3.0*xc)*xd + L*(-(rel_p*xd*(v[i,XI]-(2+cx)*xc)) + v[i,PXI]*(sx*xc + τx*v[i,PXI]/(rel_p*ωx^2))) -
                v[i,PXI]*(v[i,PXI]*sx + 2*(cx-1)*rel_p*xc)*τx/(rel_p*ωx^2)
        a = -k0*v[i,PYI]*sx2*ξ*(g*v[i,PXI]*sx*ωx + rel_p*(2+g*(v[i,XI]*(1+cx)+xc*(1-cx)))*ωx)*
            (cx2*v[i,PXI] + rel_p*sx2*xd*τx*ωx)/(2*rel_p^3*ωx^2*pry)
        b = -k0*L*υ/(2*rel_p^2*pry)
        b = b + (-g*k0*υ*(sx*xd + L*xc + (cx-1)*v[i,PXI]*τx/(rel_p*ωx^2))/(2*rel_p^2*pry))
        b = b + (-k0*υ*(cx^2*v[i,PXI]*rel_p*xd + v[i,PXI]*rel_p*xd*τx*(-τx + sx^2*ωx^2) + 
                 L*(v[i,PXI]^2 - rel_p^2*xd^2*τx*ωx^2) + cx*sx*(v[i,PXI]^2 + rel_p^2*xd^2*τx*ωx^2))/(8.0*rel_p^2*pry^3))
        b = b + (sigma*g*k0^2*v[i,PYI]^2*ξ^2/(4.0*rel_p^4*pry))
        b = b - 0.5*g*L
        cc = -k0*v[i,PYI]*ξ*(L + g*sx*xd + g*L*xc + (τx*(cx-1)*g*v[i,PXI]/(rel_p*ωx^2)))/(2*rel_p^2)
        cc = cc + (-sigma*g*k0^2*v[i,PYI]*ξ*υ/(4.0*rel_p^4*pry^2))
    end
    
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)

    # Update coordinates
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices!(K0, L, gamma_0, g, e1=nothing, e2=nothing)
    LinearTracking.linear_coast_uncoupled!(i, b::BunchView, mx, my, r56, d, t)
end

# Quadrupole
@makekernel fastgtpsa=true function magnus_quadrupole!(i, b::BunchView, K1, γ, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    k1_norm = K1 / rel_p 
    ω = sqrt(abs(k1_norm))
    s = sin(ω*L)
    sh = sinh(ω*L)
    c = cos(ω*L)
    ch = cosh(ω*L)

    if k1_norm > 0
        cx = c
        sx = s/ω
        cy = ch
        sy = sh/ω
        τy = 1
    else
        cy = c
        sy = s/ω
        cx = ch
        sx = sh/ω
        τy = -1
    end

    a = 0.5*χ*k1_norm*(sy*v[i,YI] + τy*(cy-1)/(ω^2)*v[i,PYI]/rel_p)
    b = 0.5*χ*k1_norm*(sx*v[i,XI] + τy*(1-cx)/(ω^2)*v[i,PXI]/rel_p)
    cc = v[i,PXI]*sx*(v[i,PYI]*sy + cy*rel_p*v[i,YI]) + rel_p*v[i,XI]*(cx*v[i,PYI]*sy + rel_p*(cx*cy-1)*v[i,YI])
    cc = cc * (-0.5*ξ*k1_norm/rel_p^2)
    cc = cc + 0.25*χ^2/rel_p^2*((cx-cy)*v[i,PXI]*v[i,PYI] + τy*ω^2*
         (v[i,PXI]*v[i,PYI]*sx*sy + v[i,PXI]*rel_p*(cy*sx-sy)*v[i,YI] + rel_p*v[i,XI]*(v[i,PYI]*(cx*sy-sx) + (cx*cy-1)*rel_p*v[i,YI])))
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)

    # Update coordinates
    ExactTracking.quadrupole_matrix!(i, b::BunchView, K1, L)
end

# Sextupole
@makekernel fastgtpsa=true function magnus_thick_sextupole!(i, b::BunchView, K2N, K2S, γ, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    mm = @SArray [3]
    K2N = @SArray [K2N * L/2]
    K2S = @SArray [K2S * L/2]

    # First kick of a kick-drift-kick split
    ExactTracking.multipole_kick!(i, b, mm, K2N, K2S)

    # Magnus spin rotation
    pl = sqrt(rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2)
    a = v[i,PXI]*ξ*(L^2*(-3.0*v[i,PXI]^2*v[i,PYI] + v[i,PYI]^3) - 
        3.0*L*pl*(2*v[i,PXI]*v[i,PYI]*v[i,XI] + v[i,YI]*(v[i,PXI]^2 - v[i,PYI]^2)) -
        3.0*pl^2*(2*v[i,PXI]*v[i,XI]*v[i,YI] + v[i,PYI]*(v[i,XI]^2 - v[i,YI]^2)))
    a = a + rel_p^2*χ*(L*v[i,PYI]*(2*L*v[i,PXI] + 3.0*pl*v[i,XI]) + 3.0*pl*(L*v[i,PXI] + 2*pl*v[i,XI])*v[i,YI])
    a = a * k2*L/(12.0*pl^3*rel_p^2)
    b = v[i,PYI]*ξ*(L^2*(-3.0*v[i,PXI]^2*v[i,PYI] + v[i,PYI]^3) -
        3.0*L*pl*(2*v[i,PXI]*v[i,PYI]*v[i,XI] + v[i,YI]*(v[i,PXI]^2 - v[i,PYI]^2)) -
        3.0*pl^2*(2*v[i,PXI]*v[i,XI]*v[i,YI] + v[i,PYI]*(v[i,XI]^2 - v[i,YI]^2)))
    b = b + rel_p^2*χ*(L^2*(v[i,PXI]^2 - v[i,PYI]^2) + 3.0*pl^2*(v[i,XI]^2 - v[i,YI]^2) + 3.0*L*pl*(v[i,PXI]*v[i,XI] - v[i,PYI]*v[i,YI]))
    b = b * k2*L/(12.0*pl^3*rel_p^2)
    cc = L^2*(-3.0*v[i,PXI]^2*v[i,PYI] + v[i,PYI]^3) - 3.0*L*pl*(2*v[i,PXI]*v[i,PYI]*v[i,XI] + v[i,YI]*(v[i,PXI]^2 - v[i,PYI]^2)) -
         3.0*pl^2*(2*v[i,PXI]*v[i,XI]*v[i,YI] + v[i,PYI]*(v[i,XI]^2 - v[i,YI]^2))
    cc = cc * k2*L*ξ/(12.0*pl^2*rel_p^2)
    ζ = sqrt(a^2 + b^2 + cc^2)
    sc = sincu(ζ)
    b.q *= Quaternion(-cos(ζ), a*sc, b*sc, cc*sc)

    # Drift and second kick
    ExactTracking.exact_drift!(   i, b, beta_0, gamsqr_0, tilde_m, L)
    ExactTracking.multipole_kick!(i, b, mm, K2N, K2S)
end

end