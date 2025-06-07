"""
    FieldTracking

Module implementing particle tracking through arbitrary electromagnetic fields using DifferentialEquations.jl.
"""

# Define the Field tracking method
struct Field end

# Number of temporaries needed for a single particle
MAX_TEMPS(::Field) = 0

module FieldTracking
using ..BeamTracking
using ..BeamTracking: @makekernel
using SciMLBase
const TRACKING_METHOD = Field

"""
    field_system!(du, u, p, t)

Define the ODE system for particle motion in an electromagnetic field.

# Arguments
- `du`: Vector of derivatives
- `u`: State vector [x, px, y, py, z, pz]
- `p`: Parameters tuple containing (field_func, params)
- `field_func`: Function that returns the field. Must be of the form `field_func(u, t, params)`.
    Return value should be [px, Ex, py, Ey, pz, Ez].
- `t`: Time variable
"""
function field_system!(du, u, p, t)
    field_func, params = p
    du .= field_func(u, t, params)
end

"""
    field_track!(i, v, work, L, field_func, field_params, solver, solver_params)

Track a particle through a drift space with arbitrary field using DifferentialEquations.jl.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `L`: Drift length
- `field_func`: Function that returns the field at a given position (x, y, z)
- `field_params`: Additional parameters for the field function
- `solver`: ODE solver to use
- `solver_params`: Additional parameters for the solver
"""
@makekernel function field_track!(i, v, work, L, field_func, field_params, solver, solver_params)
    @inbounds begin
        # Initial state vector
        u0 = view(v, i, :)
        
        # Set up and solve the ODE
        prob = ODEProblem(field_system!, u0, (0.0, L), (field_func, field_params))
        sol = solve(prob, solver; reltol=1e-8, abstol=1e-8, solver_params...)
        
        # Update final coordinates by assigning each component
        u0 .= sol.u[end]
    end
    return v
end

end