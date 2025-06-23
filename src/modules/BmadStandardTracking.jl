struct BmadStandard end

module BmadStandardTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView, quat_mult!
const TRACKING_METHOD = BmadStandard


# Drift - no spin rotation
@makekernel fastgtpsa=true function magnus_drift!(i, b::BunchView, β0, gamsqr_0, tilde_m, L)
    ExactTracking.exact_drift!(i, b::BunchView, β0, gamsqr_0, tilde_m, L)
end

# Solenoid -- second order Magnus expansion
@makekernel fastgtpsa=true function magnus_solenoid!(i, b::BunchView, Ks, β0, gamsqr_0, tilde_m, G, L)
    v = b.v

    if !isnothing(b.q)
        # Constants
        rel_p = 1 + v[i,PZI]
        γ = sqrt(1 + (rel_p / tilde_m)^2)
        χ = 1 + G * γ
        ξ = G * (γ - 1)
        ι = 2 * γ - 1

        # Compute effective longitudinal momentum
        pr = sqrt(rel_p^2 - (v[i,PXI] + v[i,YI]*Ks/2)^2 - (v[i,PYI] - v[i,XI]*Ks/2)^2)
        θ = Ks * L / pr
        c = cos(θ)
        s = sin(θ)

        # --- Spin rotation coefficients ---
        # A: x-plane spin rotation
        a1 = 2 * pr * (c - 1) * v[i,PYI]
        a2 = -2 * pr * s * v[i,PXI]
        a3 = Ks^2 * L * v[i,YI]
        a4 = -Ks * (2 * L * v[i,PXI] + (c - 1) * pr * v[i,XI] + pr * s * v[i,YI])
        A = (a1 + a2 + a3 + a4) * ξ / (8.0 * rel_p^2)

        # B: y-plane spin rotation
        b1 = Ks * L * (2 * v[i,PYI] + Ks * v[i,XI])
        b2 = pr * (2 * (c - 1) * v[i,PXI] + 2 * s * v[i,PYI] + c * Ks * v[i,YI] - Ks * (s * v[i,XI] + v[i,YI]))
        B = -(b1 + b2) * ξ / (8.0 * rel_p^2)

        # CC: z-plane spin rotation
        cc1 = 4.0 * pr * (v[i,PXI]^2 + v[i,PYI]^2) * s * (1 + G * ι)
        cc2 = Ks^3 * L * (v[i,XI]^2 + v[i,YI]^2) * (1 + G * ι)
        cc3 = -Ks^2 * pr * s * (v[i,XI]^2 + v[i,YI]^2) * (1 + G * ι)
        cc4 = 4.0 * Ks * ((c - 1) * pr * (v[i,PXI] * v[i,XI] + v[i,PYI] * v[i,YI]) * (1 + G * ι) +
                        L * (v[i,PXI]^2 + v[i,PYI]^2 + 4.0 * rel_p^2 + G * (4.0 * rel_p^2 + (v[i,PXI]^2 + v[i,PYI]^2) * ι)))
        cc5 = 2 * Ks^2 * L * (
            v[i,XI] * (-2 * v[i,PYI] + 2 * v[i,PXI] * s + Ks * v[i,XI]) +
            2 * (v[i,PXI] + v[i,PYI] * s) * v[i,YI] +
            Ks * v[i,YI]^2
        )
        cc6 = pr * (
            4.0 * s * (v[i,PXI]^2 + v[i,PYI]^2) +
            4.0 * Ks * v[i,PYI] * (s * v[i,XI] - 2 * v[i,YI]) -
            4.0 * Ks * v[i,PXI] * (2 * v[i,XI] + s * v[i,YI]) -
            3.0 * Ks^2 * s * (v[i,XI]^2 + v[i,YI]^2)
        )
        cc7 = c * Ks * (
            8.0 * pr * (v[i,PXI] * v[i,XI] + v[i,PYI] * v[i,YI]) +
            L * (-4.0 * (v[i,PXI]^2 + v[i,PYI]^2) + Ks^2 * (v[i,XI]^2 + v[i,YI]^2))
        )
        CC  = (cc1 + cc2 + cc3 + cc4) / (32.0 * rel_p^3)
        CC -= (cc5 + cc6 + cc7) * pr * ξ^2 / (64.0 * rel_p^4)

        # Quaternion update for spin
        ζ = sqrt(A^2 + B^2 + CC^2)
        sc = sincu(ζ)
        quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
    end

    # Update coordinates
    ExactTracking.exact_solenoid!(i, b::BunchView, Ks, β0, gamsqr_0, tilde_m, L)
end

# SBend
@makekernel fastgtpsa=true function magnus_sbend!(i, b::BunchView, g, K0, γ0, beta_gamma_0, G, L)
    
    if !isnothing(b.q)
        v = b.v
        rel_p = 1 + v[i,PZI]
        γ = sqrt(1 + (beta_gamma_0 * rel_p)^2)
        χ = 1 + G * γ
        ξ = G * (γ - 1)
        kx = g * K0

        ωx = sqrt(abs(kx)/rel_p)
        xc = abs(kx) > 0 ? (g*rel_p - K0)/kx : 0.0
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
        # Compute auxiliary variables for spin rotation
        υ = v[i, PYI]^2 * ξ - rel_p^2 * χ
        pry = sqrt(rel_p^2 - v[i, PYI]^2)

        if g == 0
            # Straight dipole
            A = -0.5 * K0 * v[i, PYI] * ξ * (L * v[i, PXI] - 0.5 * K0 * L^2) / (rel_p^2 * pry)
            B = -0.5 * K0 * (L * v[i, PXI]^2 - K0 * L^2 * v[i, PXI] + K0^2 * L^3 / 3.0) * υ / (2 * rel_p^2 * pry^3) -
                0.5 * K0 * L * υ / (rel_p^2 * pry)
            CC = -0.5 * K0 * v[i, PYI] * ξ * L / (rel_p^2)
        elseif K0 == 0
            # Bend no field
            A = 0.0
            B = -0.5 * g * L
            CC = 0.0
        else
            # Bend with field
            σ = rel_p * sx * (v[i, XI] - 3.0 * xc) * xd +
                    L * (-(rel_p * xd * (v[i, XI] - (2 + cx) * xc)) +
                    v[i, PXI] * (sx * xc + τx * v[i, PXI] / (rel_p * ωx^2))) -
                    v[i, PXI] * (v[i, PXI] * sx + 2 * (cx - 1) * rel_p * xc) * τx / (rel_p * ωx^2)

            # A: x-plane spin rotation
            A = -K0 * v[i, PYI] * sx2 * ξ *
                (g * v[i, PXI] * sx * ωx + rel_p * (2 + g * (v[i, XI] * (1 + cx) + xc * (1 - cx))) * ωx) *
                (cx2 * v[i, PXI] + rel_p * sx2 * xd * τx * ωx) /
                (2 * rel_p^3 * ωx^2 * pry)

            # B: y-plane spin rotation
            B = -0.5 * K0 * L * υ / (rel_p^2 * pry)
            B += -g * K0 * υ * (sx * xd + L * xc + (cx - 1) * v[i, PXI] * τx / (rel_p * ωx^2)) / (2 * rel_p^2 * pry)
            B += -K0 * υ * (
                    cx^2 * v[i, PXI] * rel_p * xd +
                    v[i, PXI] * rel_p * xd * τx * (-τx + sx^2 * ωx^2) +
                    L * (v[i, PXI]^2 - rel_p^2 * xd^2 * τx * ωx^2) +
                    cx * sx * (v[i, PXI]^2 + rel_p^2 * xd^2 * τx * ωx^2)
                ) / (8.0 * rel_p^2 * pry^3)
            B += σ * g * K0^2 * v[i, PYI]^2 * ξ^2 / (4.0 * rel_p^4 * pry)
            B += -0.5 * g * L

            # cc: z-plane spin rotation
            CC = -K0 * v[i, PYI] * ξ * (L + g * sx * xd + g * L * xc + (τx * (cx - 1) * g * v[i, PXI] / (rel_p * ωx^2))) / (2 * rel_p^2)
            CC += -σ * g * K0^2 * v[i, PYI] * ξ * υ / (4.0 * rel_p^4 * pry^2)
        end

        ζ = sqrt(A^2 + B^2 + CC^2)
        sc = sincu(ζ)
        quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
    end

    # Update coordinates
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, γ0)
    LinearTracking.linear_coast_uncoupled!(i, b::BunchView, mx, my, r56, d, t)
end

# Quadrupole
@makekernel fastgtpsa=true function magnus_quadrupole!(i, b::BunchView, K1, beta_gamma_0, gamsqr_0, gamma_0, beta_0, G, L)
    if !isnothing(b.q)
        v = b.v
        rel_p = 1 + v[i,PZI]
        γ = sqrt(1 + (beta_gamma_0 * rel_p)^2)
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

        # Compute spin rotation coefficients for quadrupole
        # A: x-plane spin rotation
        A = 0.5 * χ * k1_norm * (
            sy * v[i, YI] +
            τy * (cy - 1) / (ω^2) * v[i, PYI] / rel_p
        )

        # B: y-plane spin rotation
        B = 0.5 * χ * k1_norm * (
            sx * v[i, XI] +
            τy * (1 - cx) / (ω^2) * v[i, PXI] / rel_p
        )

        # CC: z-plane spin rotation
        CC = v[i, PXI] * sx * (v[i, PYI] * sy + cy * rel_p * v[i, YI]) +
            rel_p * v[i, XI] * (cx * v[i, PYI] * sy + rel_p * (cx * cy - 1) * v[i, YI])
        CC *= (-0.5 * ξ * k1_norm / rel_p^2)
        CC += 0.25 * χ^2 / rel_p^2 * (
            (cx - cy) * v[i, PXI] * v[i, PYI] +
            τy * ω^2 * (
                v[i, PXI] * v[i, PYI] * sx * sy +
                v[i, PXI] * rel_p * (cy * sx - sy) * v[i, YI] +
                rel_p * v[i, XI] * (
                    v[i, PYI] * (cx * sy - sx) +
                    (cx * cy - 1) * rel_p * v[i, YI]
                )
            )
        )

        # Quaternion update for spin
        ζ = sqrt(A^2 + B^2 + CC^2)
        sc = sincu(ζ)
        quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
    end

    # Update coordinates
    ExactTracking.quadrupole_matrix!(i, b::BunchView, K1, beta_gamma_0, gamsqr_0, gamma_0, beta_0, L)
end

# Sextupole
@makekernel fastgtpsa=true function magnus_thick_sextupole!(i, b::BunchView, K2, β0, gamsqr_0, tilde_m, G, L)
    mm = [3]
    K2N = [K2 * L/2]
    K2S = [0]

    # First kick of a kick-drift-kick split
    ExactTracking.multipole_kick!(i, b, mm, K2N, K2S)

    # ========== Magnus spin rotation ==========
    if !isnothing(b.q)
        v = b.v
        rel_p = 1 + v[i,PZI]
        γ = sqrt(1 + (rel_p / tilde_m)^2)
        χ = 1 + G * γ
        ξ = G * (γ - 1)
        pl = sqrt(rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2)

        # Common terms for reuse
        L2 = L^2
        pl2 = pl^2
        pl3 = pl^3
        rel_p2 = rel_p^2

        # Convenience pointers
        Px = v[i,PXI]
        Py = v[i,PYI]
        X  = v[i,XI]
        Y  = v[i,YI]

        # A
        term_a1 = L2 * (-3Px^2 * Py + Py^3)
        term_a2 = 3L * pl * (2Px * Py * X + Y * (Px^2 - Py^2))
        term_a3 = 3pl2 * (2Px * X * Y + Py * (X^2 - Y^2))
        term_a4 = rel_p2 * χ * (L * Py * (2L * Px + 3pl * X) + 3pl * (L * Px + 2pl * X) * Y)

        A = Px * ξ * (term_a1 - term_a2 - term_a3) + term_a4
        A *= K2 * L / (12.0 * pl3 * rel_p2)

        # B
        term_b1 = L2 * (-3Px^2 * Py + Py^3)
        term_b2 = 3L * pl * (2Px * Py * X + Y * (Px^2 - Py^2))
        term_b3 = 3pl2 * (2Px * X * Y + Py * (X^2 - Y^2))
        term_b4 = rel_p2 * χ * (L2 * (Px^2 - Py^2) + 3pl2 * (X^2 - Y^2) + 3L * pl * (Px * X - Py * Y))

        B = Py * ξ * (term_b1 - term_b2 - term_b3) + term_b4
        B *= K2 * L / (12.0 * pl3 * rel_p2)

        # CC
        term_c1 = L2 * (-3Px^2 * Py + Py^3)
        term_c2 = 3L * pl * (2Px * Py * X + Y * (Px^2 - Py^2))
        term_c3 = 3pl2 * (2Px * X * Y + Py * (X^2 - Y^2))

        CC = (term_c1 - term_c2 - term_c3) * K2 * L * ξ / (12.0 * pl2 * rel_p2)

        # Quaternion update
        ζ = sqrt(A^2 + B^2 + CC^2)
        sc = sincu(ζ)
        quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
    end
    # ========== End of Magnus spin rotation ==========

    # Drift and second kick of a kick-drift-kick split
    ExactTracking.exact_drift!(   i, b, β0, gamsqr_0, tilde_m, L)
    ExactTracking.multipole_kick!(i, b, mm, K2N, K2S)
end


# Octupole
@makekernel fastgtpsa=true function thick_octupole!(i, b::BunchView, K3N, K3S, β0, gamsqr_0, tilde_m, L)
    mm = [4]
    K3N = [K3N * L/2]
    K3S = [K3S * L/2]

    # Kick-drift-kick split
    ExactTracking.multipole_kick!(i, b, mm, K3N, K3S)
    ExactTracking.exact_drift!(   i, b, β0, gamsqr_0, tilde_m, L)
    ExactTracking.multipole_kick!(i, b, mm, K3N, K3S)
end

end