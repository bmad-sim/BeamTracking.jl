using StaticArrays

@inline function quat_mul!(a, b)
    """
    Computes the Hamilton product of a and b.
    """
    k = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
    # k[1] = a[1]*b[1] - a[2]*b[2] - a[3]*b[3] - a[4]*b[4]
    # k[2] = a[1]*b[2] + a[2]*b[1] + a[3]*b[4] - a[4]*b[3]
    # k[3] = a[1]*b[3] - a[2]*b[4] + a[3]*b[1] + a[4]*b[2]
    # k[4] = a[1]*b[4] + a[2]*b[3] - a[3]*b[2] + a[4]*b[1]

    k[1] = (a[1] * b[1] - a[3] * b[3]) - (a[2] * b[2] + a[4] * b[4])
    k[2] = (a[1] * b[2] + a[2] * b[1]) + (a[3] * b[4] - a[4] * b[3])
    k[3] = (a[1] * b[3] + a[3] * b[1]) + (a[4] * b[2] - a[2] * b[4])
    k[4] = (a[1] * b[4] + a[4] * b[1]) + (a[2] * b[3] - a[3] * b[2])

    b[1] = k[1]
    b[2] = k[2]
    b[3] = k[3]
    b[4] = k[4]
    return
end

@inline function rotate_quaternion!(q, v)
    """
    Rotates v by q
    """
    quat_mul!(q, v)
    q[2] = -q[2]
    q[3] = -q[3]
    q[4] = -q[4]
    quat_mul!(v, q)
    v[1] = q[1]
    v[2] = q[2]
    v[3] = q[3]
    v[4] = q[4]
    return
end
# sinc/sincu is zero when the real part is Inf and imag is finite
isinf_real(x::Real) = isinf(x)
isinf_real(x::Number) = isinf(real(x)) && isfinite(imag(x))

# sinhc/sinhcu is zero when the imag part is Inf and real is finite
isinf_imag(x::Real) = false
isinf_imag(x::Number) = isfinite(real(x)) && isinf(imag(x))

# sincu copied from Boost library and correct limit behavior added
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinc.hpp
"""
    sincu(x)

Compute the unnormalized sinc function ``\\operatorname{sincu}(x) = \\sin(x) / (x)`` 
with accuracy near the origin.
"""
sincu(x) = _sinc(float(x))
function _sinc(x::Union{T,Complex{T}}) where {T}
    if isinf_real(x)
        return zero(x)
    end

    nrm = Base.Math.fastabs(x)
    if nrm >= 3.3 * sqrt(sqrt(eps(T)))
        return sin(x) / x
    else
        # |x| < (eps*120)^(1/4)
        return 1 - x * x / 6
    end
end

@inline function S_dynamics!(Y_n1, Y_n2, Y_n3, RS, dt, nu0, sigma0)
    """
    Solves for the quaterion to multiply with the current spin quaternion to rotate to next step.

    Computes exp(-i/2 x /cdot /sigma) where /sigma is the vector of Pauli matrices.
    Copied from Joseph Devlin
    """

    x = MVector{3,Float64}(0.0, 0.0, 0.0)

    x[1] = 0.5 * dt * nu0
    x[2] = (-5 / 18 * dt * sigma0 * Y_n1[2]) + (-4 / 9 * dt * sigma0 * Y_n2[2]) + (-5 / 18 * dt * sigma0 * Y_n3[2])
    x[3] = (-5 / 18 * dt * sigma0 * Y_n1[1]) + (-4 / 9 * dt * sigma0 * Y_n2[1]) + (-5 / 18 * dt * sigma0 * Y_n3[1])

    exponentiate_quaternion!(x, RS)
end

@inline function exponentiate_quaternion!(x, RS)
    n = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    c = cos(n)
    s = sincu(n)

    RS[1] = -c
    RS[2] = s * x[1]
    RS[3] = s * x[2]
    RS[4] = s * x[3]

    return
end

@inline function step_orbit!(Y0, Y, dt, beta, epsilon)
    c = cos(beta * dt)
    s = sin(beta * dt)

    Y[1] = exp(-1 * epsilon * dt) * (c * Y0[1] - s * Y0[2]) + randn() * sqrt(epsilon) * sqrt(dt)
    Y[2] = exp(-1 * epsilon * dt) * (s * Y0[1] + c * Y0[2]) + randn() * sqrt(epsilon) * sqrt(dt)
end

@inline function step!(Y, R, dt, beta, nu0, epsilon, sigma0)
    """
    Solves one step of the particle dynamics using a Magnus expansion.
    """

    Y_n1 = MVector{2,Float64}(0.0, 0.0)
    Y_n2 = MVector{2,Float64}(0.0, 0.0)
    Y_n3 = MVector{2,Float64}(0.0, 0.0)
    RS = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)

    # first GL node Y
    step_orbit!(Y, Y_n1, dt * (0.5 - 0.3872983346207417), beta, epsilon)
    # second GL node Y
    step_orbit!(Y, Y_n2, dt * (0.5), beta, epsilon)
    # third GL node Y
    step_orbit!(Y, Y_n3, dt * (0.5 + 0.3872983346207417), beta, epsilon)

    # spin solve
    S_dynamics!(Y_n1, Y_n2, Y_n3, RS, dt, nu0, sigma0)

    # s_new = qpq^{-1}
    # println("R pre: ", R)
    # R_pre = R
    # println("RS pre ", RS)
    # rotate_quaternion!(RS, R)
    # println("RS post ", RS)
    quat_mul!(RS, R)
    # println("R post: ", R)
    # println("Step Difference: ", R .- R_pre)


    Y_n1[1] = Y[1]
    Y_n1[2] = Y[2]
    step_orbit!(Y_n1, Y, dt, beta, epsilon)

    return
end

@inline function magnus_solve!(Y, R, steps, dt, beta, nu0, epsilon, sigma0, true_theta)
    """
    Solves the particle dynamics from time t0 to tf with initial conditions Y0 and R0.
    """
    for _ = 1:steps
        nu = nu0 + (0.01 * sin(nu0 * (0.05 + 0.005 * randn()) * true_theta))  # THE FILES WITHOUT '0.01sine' WERE RUN WITH 0.005 * sin((0.03 + 0.001 * randn())
        step!(Y, R, dt, beta, nu, epsilon, sigma0)
        true_theta = true_theta + dt
    end
    return
end

function single_resonance_track_magnus_6_desmond(Y, R, nu0, sigma0, beta, epsilon, steps, dt, true_theta)
    """
    Function to track the spin of a particle using a single resonsance model and a Magnus expansion solver.
    """
    thread = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if thread <= size(Y)[1]
        @inbounds begin
            thread_Y = MVector{2,Float64}(Y[thread, 1], Y[thread, 2])
            thread_R = MVector{4,Float64}(R[thread, 1], R[thread, 2], R[thread, 3], R[thread, 4])

            magnus_solve!(thread_Y, thread_R, steps, dt, beta, nu0, epsilon, sigma0, true_theta)

            Y[thread, 1] = thread_Y[1]
            Y[thread, 2] = thread_Y[2]

            R[thread, 1] = thread_R[1]
            R[thread, 2] = thread_R[2]
            R[thread, 3] = thread_R[3]
            R[thread, 4] = thread_R[4]
        end
    end
    return
end


function single_resonance_track_magnus_6_threads(Y, R, nu0, sigma0, beta, epsilon, steps, dt, true_theta)
    """
    Function to track the spin of a particle using a single resonsance model and a Magnus expansion solver.
    """
    n = size(Y, 1)
    Threads.@threads for thread in 1:n
        @inbounds begin
            thread_Y = MVector{2,Float64}(Y[thread, 1], Y[thread, 2])
            thread_R = MVector{4,Float64}(R[thread, 1], R[thread, 2], R[thread, 3], R[thread, 4])

            magnus_solve!(thread_Y, thread_R, steps, dt, beta, nu0, epsilon, sigma0)

            Y[thread, 1] = thread_Y[1]
            Y[thread, 2] = thread_Y[2]

            R[thread, 1] = thread_R[1]
            R[thread, 2] = thread_R[2]
            R[thread, 3] = thread_R[3]
            R[thread, 4] = thread_R[4]
        end
    end
    return
end