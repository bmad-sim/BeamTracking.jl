"""
    RungeKuttaFieldTracking

Module implementing particle tracking through arbitrary electromagnetic fields using a 4th order Runge-Kutta method.
"""

# Define the RungeKutta tracking method
struct RungeKutta end

# Number of temporaries needed for a single particle
MAX_TEMPS(::RungeKutta) = 24  # Number of RK4 stages

module RungeKuttaTracking
using ..BeamTracking
using ..BeamTracking: @makekernel

const TRACKING_METHOD = RungeKutta

"""
    rk4_step!(u, h, field_func, params, work, i)

Perform a single 4th order Runge-Kutta step.

# Arguments
- `i`: Particle index
- `u`: State vector [x, px, y, py, z, pz]
- `work`: Work matrix (n_particles × 24)
- `t`: Current time
- `h`: Step size
- `field_func`: Function that returns the field. Must be of the form `field_func(u, t, params)`.
    Return value should be [px, Ex, py, Ey, pz, Ez].
- `params`: Additional parameters for the field function
"""
function rk4_step!(i, u, work, t, h, field_func, params)
    # Get views into work matrix for RK4 stages
    k1 = view(work, i, 1:6)
    k2 = view(work, i, 7:12)
    k3 = view(work, i, 13:18)
    k4 = view(work, i, 19:24)
    
    # Intermediate stages
    k1 .= field_func(u, 0.0, params)

    k2 .= u .+ (h/2) .* k1
    k2 .= field_func(k2, h/2, params)

    k3 .= u .+ (h/2) .* k2
    k3 .= field_func(k3, h/2, params)
    
    k4 .= u .+ h .* k3
    k4 .= field_func(k4, h, params)
    
    # Final update
    u .+= (h/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    rk4_track!(i, v, work, L, field_func, params, n_steps)

Track a particle through a drift space with arbitrary field using 4th order Runge-Kutta.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix (n_particles × 24)
- `t_span`: Time span [t_start, t_end]
- `field_func`: Function that returns the field. Must be of the form `field_func(u, t, params)`.
    Return value should be [px, Ex, py, Ey, pz, Ez].
- `params`: Additional parameters for the field function
- `n_steps`: Number of integration steps
"""
@makekernel function rk4_track!(i, v, work, t_span, field_func, params, n_steps)
    @inbounds begin
        # Create a view of the particle coordinates
        u = view(v, i, :)
        
        # Integration step size
        h = (t_span[2] - t_span[1]) / n_steps

        t = t_span[1]
        # Perform integration steps
        for _ in 1:n_steps
            rk4_step!(i, u, work, t, h, field_func, params)
            t += h
        end
    end
    return v
end

end 