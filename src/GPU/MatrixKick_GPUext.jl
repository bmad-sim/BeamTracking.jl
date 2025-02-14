module MatrixKick_GPUext
using StructArrays
using BeamTracking
using BeamTracking: get_work
using CUDA
export track!

Base.@kwdef struct Drift{T}
  L::T  # drift length / m
end

Base.@kwdef struct Quadrupole{T}
  L::T    # quadrupole length / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end

# GPU Kernel for Drift Tracking
function track_drift_kernel!(x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
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

function track!(bunch::Bunch, ele::MatrixKick_GPUext.Drift)
    L = ele.L
    tilde_m    = 1 / bunch.beta_gamma_ref
    gamsqr_ref = 1 + bunch.beta_gamma_ref^2
    beta_ref   = bunch.beta_gamma_ref / sqrt(gamsqr_ref)

    v = bunch.v  # Assuming v is a StructArray with CuArrays

    # Launch GPU Kernel
    CUDA.@sync @cuda threads=256 blocks=cld(length(v.x), 256) track_drift_kernel!(
        v.x, v.y, v.z, v.px, v.py, v.pz, L, beta_ref, tilde_m, gamsqr_ref
    )

    return bunch
end

# GPU Kernel for Quadrupole Matrix Kick
function track_quad_mx_kernel!(x, y, z, px, py, pz, k2_num, s)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
        p = 1 + pz[i]
        k2 = k2_num / p
        ks = sqrt(abs(k2)) * s
        xp = px[i] / p
        yp = py[i] / p

        focus   = k2_num >= 0
        defocus = k2_num < 0

        cx = focus * cos(ks) + defocus * cosh(ks)
        cy = focus * cosh(ks) + defocus * cos(ks)
        sx = focus * sinc(ks) + defocus * sinh(ks)
        sy = focus * sinh(ks) + defocus * sinc(ks)

        x[i]  = x[i]  * cx + xp * s * sx
        px[i] = px[i] * cx - x[i] * p * k2 * s * sx
        y[i]  = y[i]  * cy + yp * s * sy
        py[i] = py[i] * cy + y[i] * p * k2 * s * sy
        z[i] -= (s / 4) * (xp^2 + yp^2 - k2 * x[i]^2 + k2 * y[i]^2)
    end
    return nothing
end

# GPU Kernel for Quadrupole Kick
function track_quad_k_kernel!(x, y, z, px, py, pz, betgam_ref, s)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
        tilde_m = 1 / betgam_ref
        beta_ref = betgam_ref / sqrt(1 + betgam_ref^2)
        gamsqr_ref = 1 + betgam_ref^2

        p    = 1 + pz[i]
        ptr2 = px[i]^2 + py[i]^2
        ps   = sqrt(p^2 - ptr2)

        x[i] += s * px[i] / p * ptr2 / (ps * (p + ps))
        y[i] += s * py[i] / p * ptr2 / (ps * (p + ps))
        z[i] -= s * ( (1.0 + pz[i])
                     * (ptr2 - pz[i] * (2 + pz[i]) / gamsqr_ref)
                     / ( beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) * ps
                         * (beta_ref * sqrt((1.0 + pz[i])^2 + tilde_m^2) + ps))
                     - ptr2 / (2 * (1 + pz[i])^2)
                    )
    end
    return nothing
end

function track!(bunch::Bunch, ele::MatrixKick_GPUext.Quadrupole)
    L = ele.L
    k2_num = ele.Bn1 / brho(massof(bunch.species), bunch.beta_gamma_ref, chargeof(bunch.species))

    v = bunch.v
    @captured begin
    # First Matrix Kick
    CUDA.@sync @cuda threads=256 blocks=cld(length(v.x), 256) track_quad_mx_kernel!(
        v.x, v.y, v.z, v.px, v.py, v.pz, k2_num, L / 2
    )

    # Quadrupole Kick
    CUDA.@sync @cuda threads=256 blocks=cld(length(v.x), 256) track_quad_k_kernel!(
        v.x, v.y, v.z, v.px, v.py, v.pz, bunch.beta_gamma_ref, L
    )

    # Second Matrix Kick
    CUDA.@sync @cuda threads=256 blocks=cld(length(v.x), 256) track_quad_mx_kernel!(
        v.x, v.y, v.z, v.px, v.py, v.pz, k2_num, L / 2
    )
    end

    return bunch
end

end # module MatrixKick_GPUext

