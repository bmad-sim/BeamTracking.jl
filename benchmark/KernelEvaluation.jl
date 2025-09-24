using BeamTracking
using BeamTracking: get_N_particle, runkernels!, MAX_TEMPS, KernelCall, KernelChain, BunchView
using BenchmarkTools
using SciMLBase, OrdinaryDiffEq
using StaticArrays

"""
    evaluate_kernel_performance(bunch, kernel, args...; n_runs=10)

Evaluate the performance of any tracking kernel and return a dictionary of metrics.

# Arguments
- `bunch`: Initial particle bunch
- `kernel`: The kernel function to evaluate
- `args...`: Arguments to pass to the kernel
- `n_runs`: Number of runs for performance evaluation (default: 10)

# Returns
A dictionary containing the following metrics:
- `min_time`: Minimum tracking time (in nanoseconds)
- `min_memory`: Minimum memory allocation (in bytes)
- `min_allocs`: Minimum number of allocations
- `success`: Boolean whether the tracking was successful

"""
function evaluate_kernel_performance(bunch, kernel, args...; n_runs=10)
    try
        # Create kernel chain
        kc = (KernelCall(kernel, args),)
        
        # Benchmark the tracking with specified sample size and time budget
        result = @benchmark begin
            runkernels!(nothing, $bunch, $kc)
        end samples=n_runs seconds=10
        
        metrics = Dict(
            "min_time" => time(minimum(result)),
            "min_memory" => memory(minimum(result)),
            "min_allocs" => allocs(minimum(result)),
            "success" => true
        )

        return metrics
    catch e
        @warn "Tracking failed: $e"
        return Dict(
            "min_time" => NaN,
            "min_memory" => NaN,
            "min_allocs" => NaN,
            "success" => false
        )
    end 
end


function evaluate_field_track_performance(; n_runs=10, n_particles=1000, solver=Tsit5(), solver_params=(save_everystep=false,save_start=false,save_end=true,dense=false,calck=false))
    bunch = Bunch(n_particles)
    L = 1.0
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_kernel_performance(BunchView(bunch), FieldTracking.field_track!, L, field_func, params, solver, solver_params; n_runs=n_runs)
end

function evaluate_linear_track_performance(;n_runs=10, n_particles=1000)
    # suggest good default values for bunch, L, r56
    bunch = Bunch(n_particles)
    L = 1.0
    r56 = 1.0
    return evaluate_kernel_performance(BunchView(bunch), LinearTracking.linear_drift!, L, r56; n_runs=n_runs)
end

function evaluate_rk4_track_performance(;n_runs=10, n_particles=1000)
    bunch = Bunch(n_particles)
    t_span = (0.0, 1.0)
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_kernel_performance(BunchView(bunch), RungeKuttaTracking.rk4_track!, t_span, field_func, params, 10; n_runs=n_runs)
end

