struct RungeKutta end

"""
    RungeKuttaFieldTracking

Module implementing particle tracking through arbitrary electromagnetic fields using a 4th order Runge-Kutta method.
"""
module RungeKuttaTracking
using ..BeamTracking
using ..BeamTracking: @makekernel, Coords

const TRACKING_METHOD = RungeKutta

"""
    rk4_step!(u, t, h, field_func, params)

Perform a single 4th order Runge-Kutta step.

# Arguments
- `u`: State vector [x, px, y, py, z, pz]
- `t`: Current time
- `h`: Step size
- `field_func`: Function that returns the field. Must be of the form `field_func(u, t, params)`.
    Return value should be [px, Ex, py, Ey, pz, Ez].
- `params`: Additional parameters for the field function
"""
function rk4_step!(u, t, h, field_func, params)
  k1 = field_func(u, t, params)
  k2 = field_func(u .+ (h / 2) .* k1, t + h / 2, params)
  k3 = field_func(u .+ (h / 2) .* k2, t + h / 2, params)
  k4 = field_func(u .+ h .* k3, t + h, params)
  u .+= (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    rk4_track!(i, b, work, t_span, field_func, params, n_steps)

Track a particle through a drift space with arbitrary field using 4th order Runge-Kutta.

# Arguments
- `i`: Particle index
- `b`: Coords containing particle coordinates
- `t_span`: Time span [t_start, t_end]
- `field_func`: Function that returns the field. Must be of the form `field_func(u, t, params)`.
    Return value should be [px, Ex, py, Ey, pz, Ez].
- `params`: Additional parameters for the field function
- `n_steps`: Number of integration steps
"""
@makekernel function rk4_track!(i, b::Coords, t_span, field_func, params, n_steps)
  # Create a view of the particle coordinates
  u = view(b.v, i, :)

  # Integration step size
  h = (t_span[2] - t_span[1]) / n_steps

  t = t_span[1]
  # Perform integration steps
  for _ in 1:n_steps
    rk4_step!(u, t, h, field_func, params)
    t += h
  end
end

end