module MatrixKick_GPUext
export track_drift!
export track_quad!

# GPU Kernel for Drift Tracking
function track_drift!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    @inbounds begin
        ps = sqrt((1.0 + pz[i])^2 - (px[i]^2 + py[i]^2))  # P_s
        x[i] += px[i] * L / ps
        y[i] += py[i] * L / ps
        z[i] -= ( (1.0 + pz[i]) * L
                  * ((px[i]^2 + py[i]^2) - pz[i] * (2 + pz[i]) / gamsqr_ref)
                  / (beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) * ps
                     * (beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) + ps))
                )
    end
    return nothing
end

# GPU Kernel for Quadrupole Matrix Kick
function track_quad_mx!(i, x, y, z, px, py, pz, k2_num, s)
    @inbounds begin
        p = 1 + pz[i]
        k2 = k2_num / p
        ks = sqrt(abs(k2)) * s
        xp = px[i] / p
        yp = py[i] / p



        cx = k2_num >= 0 ?  cos(ks)  : cosh(ks)
        cy = k2_num >= 0 ?  cosh(ks) : cos(ks)
        sx = k2_num >= 0 ?  sinc(ks) : sinh(ks)
        sy = k2_num >= 0 ?  sinh(ks) : sinc(ks)


        x[i]  = x[i]  * cx + xp * s * sx
        px[i] = px[i] * cx - x[i] * p * k2 * s * sx
        y[i]  = y[i]  * cy + yp * s * sy
        py[i] = py[i] * cy + y[i] * p * k2 * s * sy
        z[i] -= (s / 4) * (xp^2 + yp^2 - k2 * x[i]^2 + k2 * y[i]^2)
    end
    return nothing
end

# GPU Kernel for Quadrupole Kick
function track_quad_k!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    @inbounds begin
        p    = 1 + pz[i]
        ptr2 = px[i]^2 + py[i]^2
        ps   = sqrt(p^2 - ptr2)

        x[i] += L * px[i] / p * ptr2 / (ps * (p + ps))
        y[i] += L * py[i] / p * ptr2 / (ps * (p + ps))
        z[i] -= L * ( (1.0 + pz[i])
                     * (ptr2 - pz[i] * (2 + pz[i]) / gamsqr_ref)
                     / ( beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) * ps
                         * (beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) + ps))
                     - ptr2 / (2 * (1 + pz[i])^2)
                    )
    end
    return nothing
end

function track_quad!(i, x, y, z, px, py, pz, beta_ref, tilde_m, gamsqr_ref, L, Bn1, brho)
    k2_num = Bn1 / brho
    track_quad_mx!(i, x, y, z, px, py, pz, k2_num, L / 2)
    # Quadrupole Kick
    track_quad_k!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    # Second Matrix Kick
    track_quad_mx!(i, x, y, z, px, py, pz, k2_num, L / 2)

    return nothing
end

end # module MatrixKick_GPUext

