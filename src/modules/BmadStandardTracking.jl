struct BmadStandard end

module BmadStandardTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView, quat_mult!
const TRACKING_METHOD = BmadStandard


@makekernel fastgtpsa=true function bmad_cavity!(i, b::BunchView, V, wave_number, φ0, β0, gamsqr_0, tilde_m, p0c, L)
  v = b.v
  z = v[i, ZI]
  # ============= Kick - Drift - Kick scheme =============
  # First kick
  rel_p = 1 + v[i, PZI]
  φ = φ0 - 2π * wave_number * z * sqrt(rel_p^2 + tilde_m^2) / (rel_p * sqrt(1 + tilde_m^2))
  dE = V * sin(φ) / 2
  dE_rel = dE / p0c
  S = sqrt(rel_p^2 + tilde_m^2)
  Δ = (2S + dE_rel) * dE_rel
  v[i, PZI] += Δ / (rel_p + sqrt(rel_p^2 + Δ^2))

  # Drift
  ExactTracking.exact_drift!(i, b::BunchView, β0, gamsqr_0, tilde_m, L)

  # Second kick
  rel_p = 1 + v[i, PZI]
  φ += 2π * wave_number * (v[i, ZI] - z) * sqrt(rel_p^2 + tilde_m^2) / (rel_p * sqrt(1 + tilde_m^2))
  dE = V * sin(φ) / 2
  dE_rel = dE / p0c
  S = sqrt(rel_p^2 + tilde_m^2)
  Δ = (2S + dE_rel) * dE_rel
  v[i, PZI] += Δ / (rel_p + sqrt(rel_p^2 + Δ^2))
end

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
@makekernel fastgtpsa=true function magnus_sbend!(i, b::BunchView, g, K0, γ0, βγ0, G, L)
    
    if !isnothing(b.q)
        v = b.v
        rel_p = 1 + v[i,PZI]
        γ = sqrt(1 + (βγ0 * rel_p)^2)
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

            # CC: z-plane spin rotation
            CC = -K0 * v[i, PYI] * ξ * (L + g * sx * xd + g * L * xc + (τx * (cx - 1) * g * v[i, PXI] / (rel_p * ωx^2))) / (2 * rel_p^2)
            CC += -σ * g * K0^2 * v[i, PYI] * ξ * υ / (4.0 * rel_p^4 * pry^2)
        end

        ζ = sqrt(A^2 + B^2 + CC^2)
        sc = sincu(ζ)
        quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
    end

    # Update coordinates
    ExactTracking.exact_bend!(i, b::BunchView, g * L, g, K0, 1/βγ0, βγ0/γ0, L)
end

# Quadrupole
@makekernel fastgtpsa=true function magnus_quadrupole!(i, b::BunchView, K1, βγ0, tilde_m, G, L)
    v = b.v
    rel_p = 1 + v[i,PZI]
    if !isnothing(b.q)
        γ = sqrt(1 + (βγ0 * rel_p)^2)
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
    ExactTracking.quadrupole_matrix!(i, b::BunchView, K1, L)
    # beta != beta_0 correction
    rel_e = sqrt(rel_p^2 + tilde_m^2)
    v[i,ZI] += L * ( tilde_m^2 * v[i,PZI] * (2 + v[i,PZI]) / ( rel_e * ( rel_p * sqrt(1 + tilde_m^2) + rel_e ) ) )
end

@makekernel fastgtpsa=true function magnus_combined_func!(i, b::BunchView, g, k0, k1, tilde_m, G, L)
  v = b.v
  rel_p  = 1 + v[i, PZI]
  inv_rel_p = 1 / rel_p
  kx  = k1 + g * k0
  xc  = (g * rel_p - k0) / kx

  ωx  = sqrt(abs(kx)) * sqrt(inv_rel_p)
  ωy  = sqrt(abs(k1)) * sqrt(inv_rel_p)

  
  arg = ωx * L
  if arg < 1e-6
    sx = (1 - sign(kx) * arg^2 / 6) * L      # s_x
    cx = 1 - sign(kx) * arg^2 / 2            # c_x
    z2 = g * L^2 / (2 * rel_p)                        # z2
  elseif kx > 0
    sx = sin(arg) / ωx                        # s_x
    cx = cos(arg)                               # c_x
    z2 = -sign(kx) * g * (1 - cx) / (rel_p * ωx^2)# z2
  else
    sx = sinh(arg) / ωx                       # s_x
    cx = cosh(arg)                              # c_x
    z2 = -sign(kx) * g * (1 - cx) / (rel_p * ωx^2)# z2
  end

  arg = ωy * L
  if arg < 1e-6
    sy = (1 + sign(k1) * arg^2 / 6) * L            # s_y
    cy = 1 + sign(k1) * arg^2 / 2                  # c_y
  elseif k1 < 0
    sy = sin(arg) / ωy                             # s_y
    cy = cos(arg)                                    # c_y
  else
    sy = sinh(arg) / ωy                            # s_y
    cy = cosh(arg)                                   # c_y
  end

  x0  = v[i,  XI] 
  px0 = v[i, PXI]
  y0  = v[i,  YI]
  py0 = v[i, PYI]


  if !isnothing(b.q)
    γ = sqrt(1 + (rel_p / tilde_m)^2)
    χ = 1 + G * γ
    ξ = G * (γ - 1)
    ν = 2 * (1 + G) - χ

    η = ωx^2 + ωy^2
    μ = ωx^2 - ωy^2

    c = cx * cy - 1
    xd = x0 - xc

    a1 = b1 = cc1 = a2 = b2 = cc2 = 0.0

    if kx > 0 && k1 > 0 
          a1 = k0 * ξ * ωy^2 * (
              py0 * ωx^2 * (-cy * px0 * sx - c * rel_p * xd * ωx) -
              (cx * px0 * py0 * sy + c * px0 * rel_p * y0 + px0 * rel_p * sx * sy * y0 * ωx^2 -
              rel_p * xd * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωx^3) * ωy^2
          )
          a1 += k1 * χ * (
              (cy - 1) * py0 * rel_p * (1 + g * xc) * η +
              g * py0 * (px0 * (cy * sx - cx * sy) + cx * cy * rel_p * xd + rel_p * xd * (sx * sy * ωx^2 - 1)) * ωy^2 +
              rel_p * y0 * ωy^2 * (
                  rel_p * sy * η +
                  g * px0 * (sx * sy * ωy^2 - c) +
                  g * rel_p * (sy * xc * η + cy * sx * xd * ωx^2 + cx * sy * xd * ωy^2)
              )
          )
          a1 *= 0.5 / (rel_p^3 * η * ωy^2)

          b1 = 2 * g * px0 * (
              k1 * L * px0 - 2 * (cx - 1) * k0 * rel_p +
              k1 * (rel_p * (x0 + 3 * xc) - cx * (px0 * sx + 4 * rel_p * xc) - cx^2 * rel_p * xd)
          ) * χ
          b1 += k0 * py0 * (L * py0 + cy * py0 * sy + (cy^2 - 1) * y0) * ν * ωx^2
          b1 -= 2 * g * rel_p * (
              -sx * xd * (2 * k0 * rel_p + k1 * px0 * sx + k1 * rel_p * (4 * xc + cx * xd)) * χ +
              L * rel_p * (2 * rel_p - (2 * k0 * xc + k1 * (x0^2 - 2 * x0 * xc + 3 * xc^2)) * χ)
          ) * ωx^2
          b1 += 4 * k1 * rel_p * χ * (px0 * (1 - cx) + rel_p * (L * xc + sx * xd) * ωx^2)
          b1 += k0 * ωx^2 * (
              L * (4 + px0^2) * χ + χ * (
                  8 * L * (rel_p - 1) + 4 * L * (rel_p - 1)^2 +
                  cx * px0^2 * sx + (cx^2 - 1) * px0 * rel_p * xd * ωx -
                  px0 * rel_p * sx^2 * xd * ωx^3 + rel_p^3 * (L - cx * sx) * xd^2 * ωx^4
              )
          )
          b1 -= rel_p^2 * (L - cy * sy) * y0^2 * ν * ωy^2
          b1 += py0 * y0 * ν * ((cy^2 - 1) * (rel_p - 1) + rel_p * sy^2 * ωy^2)
          b1 /= 8 * rel_p^3 * ωx^2

          cc1 = k0 * rel_p * (1 + g * xc) * (py0 * sy + (cy - 1) * rel_p * y0) * η
          cc1 += kx * py0 * (-c * px0 + cy * rel_p * sx * xd * ωx^2)
          cc1 += kx * (
              px0 * py0 * sx * sy + px0 * rel_p * (cy * sx - cx * sy) * y0 +
              rel_p * xd * (cx * py0 * sy + cx * cy * rel_p * y0 + rel_p * y0 * (sx * sy * ωx^2 - 1))
          ) * ωy^2
          cc1 += k1 * (
              px0 * (
                  py0 * (c + sx * sy * ωx^2) + rel_p * y0 * (cy * sx * ωx^2 + cx * sy * ωy^2)
              ) +
              rel_p * (
                  py0 * (sy * xc * η - cy * sx * xd * ωx^3 + cx * sy * xd * ωx^3) +
                  rel_p * y0 * ((cy - 1) * xc * η + xd * ωx^3 * (c - sx * sy * ωy^2))
              )
          )
          cc1 *= -0.5 * ξ / (rel_p^3 * η)

          a2 = rel_p * y0 * ωy^2 * (
              (px0 - rel_p * (L * xc - 2 * sy * xc + sx * xd) * η) * ωx^2 +
              cx * px0 * (η - cy * ωx^2) +
              px0 * (-1 + 2 * sx * sy * ωx^2) * ωy^2 +
              cx * (cy * px0 + 2 * rel_p * sy * xd * ωx^2) * ωy^2 +
              cy * (-px0 * η + rel_p * ωx^2 * (-L * xc * η + sx * xd * μ))
          )
          a2 += py0 * (
              -rel_p * xc * ωx^2 * (η - 2 * cy * η + ωx^2 + L * sy * η * ωy^2) +
              ωy^2 * (
                  -rel_p * ωx^2 * (x0 + xd - 2 * cx * cy * xd - sx * sy * xd * μ) +
                  px0 * (2 * cy * sx * ωx^2 - sy * (η + cx * μ))
              )
          )
          a2 *= 0.25 * k0 * kx * ξ * χ / (rel_p^4 * η * ωx^2 * ωy^2)

          b2 = k0 * k1 * (L - sy) * ξ * χ * (py0^2 - rel_p^2 * y0^2 * ωy^2) / (4 * rel_p^4 * ωy^2)

          cc2 = px0 * (
              py0 * ((cx - cy - 1) * η + cx * cy * ωy^2 + ωx^2 * (2 - cx * cy + 2 * sx * sy * ωy^2)) +
              rel_p * y0 * ωy^2 * (2 * cy * sx * ωx^2 - sy * (η + cx * μ))
          )
          cc2 -= rel_p * ωx^2 * (
              L * xc * η * ((cy + 1) * py0 + rel_p * sy * y0 * ωy^2) +
              py0 * (
                  -2 * sy * (xc * η + cx * xd * ωy^2) + sx * xd * (η - cy * μ)
              ) +
              rel_p * y0 * (
                  2 * xc * ωx^2 - 2 * cy * (xc * η + cx * xd * ωy^2) + ωy^2 * (2 * x0 - sx * sy * xd * μ)
              )
          )
          cc2 *= 0.25 * k1 * kx * χ^2 / (rel_p^4 * η * ωx^2 * ωy^2)
    elseif kx > 0 && k1 < 0
        if abs(μ) > 1e-14
            a1 = k0 * ξ * ωy^2 * (
                py0 * ωx^2 * (-cy * px0 * sx - c * rel_p * xd * ωx) +
                (cx * px0 * py0 * sy + c * px0 * rel_p * y0 + px0 * rel_p * sx * sy * y0 * ωx^2 -
                 rel_p * xd * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωx^3) * ωy^2
            )
            a1 += k1 * χ * (
                -((cy - 1) * py0 * rel_p * (1 + g * xc) * ωx^2) +
                (rel_p * y0 * (-c * g * px0 + rel_p * (sy + g * sy * xc + cy * g * sx * xd) * ωx^2) +
                 py0 * (-rel_p - cx * g * px0 * sy + cy * (rel_p + g * (rel_p - 1) * (xc + cx * xd) + g * (px0 * sx + xc + cx * xd)) -
                        g * rel_p * (xc + xd - sx * sy * xd * ωx^2))) * ωy^2 -
                rel_p * sy * (rel_p + g * (rel_p - 1) * (xc + cx * xd) + g * (px0 * sx + xc + cx * xd)) * y0 * ωy^4
            )
            a1 *= 0.5 / (rel_p^3 * μ * ωy^2)

            cc1 = -k0 * rel_p * (1 + g * xc) * (py0 * sy + (cy - 1) * rel_p * y0) * μ
            cc1 += kx * (c * px0 * py0 + px0 * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωy^2 +
                       rel_p * xd * (-cy * py0 * sx * ωx^2 + (cx * py0 * sy + cx * cy * rel_p * y0 +
                       rel_p * y0 * (-1 + sx * sy * ωx^2)) * ωy^2))
            cc1 -= k1 * (
                px0 * (py0 * (c + sx * sy * ωx^2) + rel_p * y0 * (cy * sx * ωx^2 - cx * sy * ωy^2)) +
                rel_p * (
                    py0 * ωx^2 * (-cy * sx * xd * ωx + sy * (xc + cx * xd * ωx)) - py0 * sy * xc * ωy^2 +
                    rel_p * y0 * ((cy - 1) * xc * μ + xd * ωx^3 * (c - sx * sy * ωy^2))
                )
            )
            cc1 *= 0.5 * ξ / (rel_p^3 * μ)

            a2 = -rel_p * y0 * ωy^2 * (
                (1 + cx) * (cy - 1) * px0 * ωx^2 + px0 * ((cx - 1) * (1 + cy) + 2 * sx * sy * ωx^2) * ωy^2 +
                rel_p * ωx^2 * ((1 + cy) * L * xc - 2 * sy * xc + sx * xd - cy * sx * xd) * ωx^2 -
                ((1 + cy) * L * xc + (1 + cy) * sx * xd - 2 * sy * (xc + cx * xd)) * ωy^2
            )
            a2 += py0 * (
                -rel_p * xc * ωx^2 * μ * (-2 + 2 * cy + L * sy * ωy^2) +
                ωy^2 * (2 * cy * (px0 * sx + cx * rel_p * xd) * ωx^2 -
                       px0 * sy * ((1 + cx) * ωx^2 + (cx - 1) * ωy^2) +
                       rel_p * xd * ωx^2 * (-2 + sx * sy * η))
            )
            a2 *= 0.25 * k0 * kx * ξ * χ / (rel_p^4 * ωx^2 * ωy^2 * μ)

            cc2 = px0 * ((1 + cx) * (cy - 1) * py0 * ωx^2 + py0 * ((cx - 1) * (1 + cy) + 2 * sx * sy * ωx^2) * ωy^2 +
                         rel_p * y0 * ωy^2 * ((2 * cy * sx - (1 + cx) * sy) * ωx^2 - (1 - cx) * sy * ωy^2))
            cc2 -= rel_p * ωx^2 * (
                py0 * (2 * sy * xc + (cy - 1) * sx * xd) * ωx^2 + py0 * ((1 + cy) * sx * xd - 2 * sy * (xc + cx * xd)) * ωy^2 +
                L * xc * μ * (-((1 + cy) * py0) + rel_p * sy * y0 * ωy^2) +
                rel_p * y0 * (2 * (cy - 1) * xc * ωx^2 +
                              (x0 + xc - 2 * cy * xc + xd - 2 * cx * cy * xd - sx * sy * xd * ωx^2) * ωy^2 -
                              sx * sy * xd * ωy^4)
            )
            cc2 *= 0.25 * k1 * kx * χ^2 / (rel_p^4 * ωx^2 * ωy^2 * μ)
        else
            error("μ ≈ 0 should not happen for k0 > 0 and k1 < 0")
        end
    elseif kx < 0 && k1 > 0
        if abs(μ) > 1e-14
            a1 = k0 * ξ * ωy^2 * (
                py0 * ωx^2 * (-cy * px0 * sx - c * rel_p * xd * ωx) +
                (cx * px0 * py0 * sy + c * px0 * rel_p * y0 - px0 * rel_p * sx * sy * y0 * ωx^2 +
                 rel_p * xd * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωx^3) * ωy^2
            )
            a1 += k1 * χ * (
                (cy - 1) * py0 * rel_p * (1 + g * xc) * ωx^2 +
                (rel_p * y0 * (c * g * px0 + rel_p * (sy + g * sy * xc + cy * g * sx * xd) * ωx^2) +
                 py0 * (rel_p - cy * (rel_p + g * px0 * sx + g * rel_p * xc + cx * g * rel_p * xd) +
                        g * (cx * px0 * sy + rel_p * (xc + xd + sx * sy * xd * ωx^2)))) * ωy^2 -
                rel_p * sy * (rel_p + g * px0 * sx + g * rel_p * xc + cx * g * rel_p * xd) * y0 * ωy^4
            )
            a1 *= 0.5 / (rel_p^3 * μ * ωy^2)

            cc1 = kx * py0 * (-px0 * c - cy * rel_p * sx * xd * ωx^2)
            cc1 += kx * (
                px0 * py0 * sx * sy + px0 * rel_p * (cy * sx - cx * sy) * y0 +
                rel_p * xd * (cx * py0 * sy + cx * cy * rel_p * y0 - rel_p * y0 * (1 + sx * sy * ωx^2))
            ) * ωy^2
            cc1 -= k0 * rel_p * (1 + g * xc) * (py0 * sy + (cy - 1) * rel_p * y0) * μ
            cc1 += k1 * (
                px0 * (py0 * (c - sx * sy * ωx^2) - rel_p * y0 * (cy * sx * ωx^2 - cx * sy * ωy^2)) +
                rel_p * (
                    py0 * ωx^2 * (cy * sx * xd * ωx - sy * (xc + cx * xd * ωx)) + py0 * sy * xc * ωy^2 +
                    rel_p * y0 * (-(cy - 1) * xc * μ + xd * ωx^3 * (c - sx * sy * ωy^2))
                )
            )
            cc1 *= 0.5 * ξ / (rel_p^3 * μ)

            a2 = rel_p * y0 * ωy^2 * (
                (1 + cx) * (cy - 1) * px0 * ωx^2 + px0 * ((cx - 1) * (1 + cy) - 2 * sx * sy * ωx^2) * ωy^2 +
                rel_p * ωx^2 * (
                    2 * sy * xc + (cy - 1) * sx * xd) * ωx^2 + ((1 + cy) * sx * xd - 2 * sy * (xc + cx * xd)) * ωy^2 -
                    (1 + cy) * L * xc * μ
            )
            a2 += py0 * (
                -rel_p * xc * ωx^2 * μ * (2 - 2 * cy - L * sy * ωy^2) +
                ωy^2 * (
                    -2 * cy * (px0 * sx + cx * rel_p * xd) * ωx^2 +
                    px0 * sy * ((1 + cx) * ωx^2 + (cx - 1) * ωy^2) +
                    rel_p * xd * ωx^2 * (2 + sx * sy * η)
                )
            )
            a2 *= 0.25 * k0 * kx * ξ * χ / (rel_p^4 * ωx^2 * ωy^2 * μ)

            cc2 = px0 * (
                (1 + cx) * (cy - 1) * py0 * ωx^2 + ((cx - 1) * (1 + cy) * py0 + (-2 * py0 * sx * sy +
                rel_p * (-2 * cy * sx + sy + cx * sy) * y0) * ωx^2) * ωy^2 + (-1 + cx) * rel_p * sy * y0 * ωy^4
            )
            cc2 += rel_p * ωx^2 * (
                py0 * (2 * sy * xc + (cy - 1) * sx * xd) * ωx^2 + py0 * ((1 + cy) * sx * xd - 2 * sy * (xc + cx * xd)) * ωy^2 -
                L * xc * μ * ((1 + cy) * py0 + rel_p * sy * y0 * ωy^2) +
                rel_p * y0 * (
                    2 * (cy - 1) * xc * ωx^2 +
                    (x0 + xc - 2 * cy * xc + xd - 2 * cx * cy * xd + sx * sy * xd * ωx^2) * ωy^2 +
                    sx * sy * xd * ωy^4
                )
            )
            cc2 *= 0.25 * k1 * kx * χ^2 / (rel_p^4 * ωx^2 * ωy^2 * μ)
        else
            error("μ ≈ 0 should not happen for k0 < 0 and k1 > 0")
        end
    else
        # kx < 0 && k1 < 0
        a1 = k0 * ξ * ωy^2 * (
            py0 * ωx^2 * (-cy * px0 * sx - c * rel_p * xd * ωx) -
            (cx * px0 * py0 * sy + c * px0 * rel_p * y0 - px0 * rel_p * sx * sy * y0 * ωx^2 +
             rel_p * xd * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωx^3) * ωy^2
        )
        a1 += k1 * χ * (
            -((cy - 1) * py0 * rel_p * (1 + g * xc) * η) +
            g * py0 * (cx * px0 * sy - cy * (px0 * sx + cx * rel_p * xd) + rel_p * xd * (1 + sx * sy * ωx^2)) * ωy^2 +
            rel_p * y0 * ωy^2 * (
                rel_p * sy * η + g * px0 * (c + sx * sy * ωy^2) +
                g * rel_p * (sy * xc * η + cy * sx * xd * ωx^2 + cx * sy * xd * ωy^2)
            )
        )
        a1 *= 0.5 / (rel_p^3 * η * ωy^2)

        b1 = 2 * g * px0 * (
            2 * (cx - 1) * k0 * rel_p - k1 * (L * px0 - cx * px0 * sx + x0 + (rel_p - 1) * x0) +
            (-3 + 4 * cx) * k1 * rel_p * xc + cx^2 * k1 * rel_p * xd
        ) * χ
        b1 += k0 * py0 * (L * py0 + cy * py0 * sy + (cy^2 - 1) * y0) * ν * ωx^2
        b1 -= 2 * g * rel_p * (
            -sx * xd * (2 * k0 * rel_p + k1 * px0 * sx + k1 * rel_p * (4 * xc + cx * xd)) * χ +
            L * rel_p * (2 * rel_p - (2 * k0 * xc + k1 * (x0^2 - 2 * x0 * xc + 3 * xc^2)) * χ)
        ) * ωx^2
        b1 += 4 * k1 * rel_p * χ * ((-1 + cx) * px0 + rel_p * (L * xc + sx * xd) * ωx^2)
        b1 += k0 * ωx^2 * (
            L * (4 + px0^2) * χ + χ * (
                8 * L * (rel_p - 1) + 4 * L * (rel_p - 1)^2 + cx * px0^2 * sx +
                (cx^2 - 1) * px0 * rel_p * xd * ωx + px0 * rel_p * sx^2 * xd * ωx^3 -
                rel_p^2 * (L - cx * sx) * xd^2 * ωx^4
            )
        )
        b1 += rel_p^2 * (L - cy * sy) * y0^2 * ν * ωy^2
        b1 += py0 * y0 * ν * ((cy^2 - 1) * (rel_p - 1) - rel_p * sy^2 * ωy^2)
        b1 /= 8 * rel_p^3 * ωx^2

        cc1 = -k0 * rel_p * (1 + g * xc) * (py0 * sy + (cy - 1) * rel_p * y0) * η
        cc1 -= kx * (
            c * px0 * py0 + px0 * (py0 * sx * sy + rel_p * (cy * sx - cx * sy) * y0) * ωy^2 +
            rel_p * xd * (cy * py0 * sx * ωx^2 + (cx * py0 * sy + cx * cy * rel_p * y0 -
            rel_p * y0 * (1 + sx * sy * ωx^2)) * ωy^2)
        )
        cc1 += k1 * (
            px0 * (py0 * (c - sx * sy * ωx^2) - rel_p * y0 * (cy * sx * ωx^2 + cx * sy * ωy^2)) +
            rel_p * (
                -py0 * sy * xc * η + py0 * (cy * sx - cx * sy) * xd * ωx^3 -
                rel_p * y0 * ((cy - 1) * xc * η + xd * ωx^3 * (c + sx * sy * ωy^2))
            )
        )
        cc1 *= 0.5 * ξ / (rel_p^3 * η)

        a2 = rel_p * y0 * ωy^2 * (
            px0 * ((1 + cx - cy) * η - cx * cy * ωx^2 + (-2 + cx * cy - 2 * sx * sy * ωx^2) * ωy^2) +
            rel_p * ωx^2 * ((1 + cy) * L * xc * η + sx * xd * (η - cy * μ) - 2 * sy * (xc * η + cx * xd * ωy^2))
        )
        a2 += py0 * (
            rel_p * xc * η * ωx^2 * (-2 + 2 * cy + L * sy * ωy^2) +
            ωy^2 * (2 * cy * px0 * sx * ωx^2 - px0 * sy * (η + cx * μ) + rel_p * xd * ωx^2 * (2 * c - sx * sy * μ))
        )
        a2 *= -0.25 * k0 * kx * ξ * χ / (rel_p^4 * η * ωx^2 * ωy^2)

        b2 = -k0 * k1 * (L - sy) * ξ * χ * (py0^2 + rel_p^2 * y0^2 * ωy^2) / (4 * rel_p^4 * ωy^2)

        cc2 = px0 * (
            py0 * ((1 + cx - cy) * η - cx * cy * ωx^2 + (-2 + cx * cy - 2 * sx * sy * ωx^2) * ωy^2) +
            rel_p * y0 * ωy^2 * (-2 * cy * sx * ωx^2 + sy * (η + cx * μ))
        )
        cc2 += rel_p * ωx^2 * (
            L * xc * η * ((1 + cy) * py0 - rel_p * sy * y0 * ωy^2) +
            py0 * (sx * xd * (η - cy * ωx^2 + cy * ωy^2) - 2 * sy * (xc * η + cx * xd * ωy^2)) -
            rel_p * y0 * (2 * (cy - 1) * xc * η + xd * ωy^2 * (2 * c - sx * sy * μ))
        )
        cc2 *= 0.25 * k1 * kx * χ^2 / (rel_p^4 * η * ωx^2 * ωy^2)
    end

    # Final sum
    A = a1 + a2
    B = b1 + b2
    CC = cc1 + cc2

    ζ = sqrt(A^2 + B^2 + CC^2)
    sc = sincu(ζ)
    quat_mult!(@SVector[-cos(ζ), A*sc, B*sc, CC*sc], b.q)
  end

  x0 -= xc

  # Update transverse
  v[i,  XI] = cx * x0 + sx * px0 * inv_rel_p + xc
  v[i, PXI] = -sign(kx) * ωx^2 * rel_p * sx * x0 + cx * px0
  v[i,  YI] = cy * y0 + sy * py0 * inv_rel_p
  v[i, PYI] = sign(k1) * ωy^2 * rel_p * sy * y0 + cy * py0

  # Longitudinal update
  rel_e = sqrt(rel_p^2 + tilde_m^2)
  v[i, ZI] += L * ( tilde_m^2 * v[i,PZI] * (2 + v[i,PZI]) / ( rel_e * ( rel_p * sqrt(1 + tilde_m^2) + rel_e ) ) ) + 
          (-g * xc * L) +
          (-g * sx) * x0 +
          z2 * px0 +
          (-sign(kx) * ωx^2 * (L - cx * sx) / 4) * x0^2 +
          (sign(kx) * ωx^2 * sx^2 / (2 * rel_p)) * x0 * px0 +
          (-(L + cx * sx) / (4 * rel_p^2)) * px0^2 +
          (sign(k1) * ωy^2 * (L - cy * sy) / 4) * y0^2 +
          (-sign(k1) * ωy^2 * sy^2 / (2 * rel_p)) * y0 * py0 +
          (-(L + cy * sy) / (4 * rel_p^2)) * py0^2
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
      pl2 = rel_p^2 - v[i,PXI]^2 - v[i,PYI]^2
      if pl2 <= 0
        b.state[i] = State.Lost
        @warn "Particle lost in sextupole (transverse momentum too large)"
      else
        pl = sqrt(pl2)

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

@makekernel fastgtpsa=true function hwang_edge!(i, b::BunchView, e, k0, k1, upstream)
  cos_e = cos(e); 
  sin_e = sin(e); 
  tan_e = sin_e / cos_e; 
  sec_e = 1 / cos_e
  gt = k0 * tan_e
  t2 = tan_e * tan_e
  gt2 = k0 * t2
  gs2 = k0 * sec_e^2
  k1_tane = k1 * tan_e

  s = 2 * upstream - 1

  v = b.v
  v1_2 = v[i, XI] * v[i, XI]
  v3_2 = v[i, YI] * v[i, YI]
  e_factor = 1 / (1 + v[i, PZI])
  fg_factor = 0

  dx  = s * (-gt2 * v[i, XI]^2 + gs2 * v[i, YI]^2) * e_factor / 2

  dpx = e_factor * (s * gt2 * (v[i, XI] * v[i, PXI] - v[i, YI] * v[i, PYI]) + 
        k1_tane * (v1_2 - v3_2) + k0 * gt * (
            0.5 * upstream * (1 + 2 * t2) * v3_2 
          + 0.25 * (s - 1) * t2 * (v1_2 + v3_2)
        )
      )

  dy  = s * gt2 * v[i, XI] * v[i, YI] * e_factor

#=
  if upstream
    dpx = (gt * g_tot * (1 + 2 * tan_e^2) * v[i, YI]^2 / 2 + gt2 * (v[i, XI] * v[i, PXI] - v[i, YI] * v[i, PYI]) + k1_tane * (v[i, XI]^2 - v[i, YI]^2)) * e_factor
    dpy = (fg_factor * v[i, YI] - gt2 * v[i, XI] * v[i, PYI] - (k0 + gt2) * v[i, PXI] * v[i, YI] - 2 * k1_tane * v[i, XI] * v[i, YI]) * e_factor
  else
    dpx = (gt2 * (v[i, YI] * v[i, PYI] - v[i, XI] * v[i, PXI]) + k1_tane * (v[i, XI]^2 - v[i, YI]^2) - gt * gt2 * (v[i, XI]^2 + v[i, YI]^2) / 2) * e_factor
    dpy = (fg_factor * v[i, YI] + gt2 * v[i, XI] * v[i, PYI] + (k0 + gt2) * v[i, PXI] * v[i, YI] + (- 2 * k1_tane + gt * gs2) * v[i, XI] * v[i, YI]) * e_factor
  end
=#

  dpy = (fg_factor * v[i, YI] - s * (gt2 * v[i, XI] * v[i, PYI] + (k0 + gt2) * v[i, PXI] * v[i, YI]) + ((1 - upstream) * gt * gs2 - 2 * k1_tane) * v[i, XI] * v[i, YI]) * e_factor


  dz = e_factor^2 * 0.5 * (v[i, YI]^2 * fg_factor +
            v[i, XI]^3 * (4.0 * k1_tane - gt * gt2) / 6.0 + 0.5 * v[i, XI]*v[i, YI]^2 * (-4.0 * k1_tane + gt * gs2) +
            s * ((v[i, XI]^2*v[i, PXI] - 2.0 * v[i, XI]*v[i, YI]*v[i, PYI]) * gt2 - v[i, PXI]*v[i, YI]^2 * gs2))


  v[i, PXI] += dpx + gt * v[i, XI]
  v[i, PYI] += dpy - gt * v[i, YI]
  v[i,  XI] += dx
  v[i,  YI] += dy
  v[i,  ZI] += dz
end

end