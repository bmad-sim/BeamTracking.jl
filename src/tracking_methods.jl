# ========== Yoshida ===========================
abstract type AbstractYoshida end
macro def_integrator_struct(name)
  quote
    struct $(esc(name)) <: AbstractYoshida
      order::Int
      num_steps::Int 
      ds_step::Float64
      radiation_damping_on::Bool
      radiation_fluctuations_on::Bool
  
      function $(esc(name))(; order::Int=4, num_steps::Int=-1, ds_step::Float64=-1.0, radiation_damping_on::Bool=false, radiation_fluctuations_on::Bool=false)
        _order = order
        _num_steps = num_steps
        _ds_step = ds_step
        if _order ∉ (2, 4, 6, 8)
          error("Symplectic integration only supports orders 2, 4, 6, and 8")
        elseif _num_steps == -1 && _ds_step == -1.0
          _num_steps = 1
        elseif _num_steps > 0 && _ds_step > 0
          error("Only one of num_steps or ds_step should be specified")
        elseif _num_steps < 1 && _ds_step <= 0
          error("Invalid step size")
        elseif _num_steps > 0
          _ds_step = -1.0
        elseif _ds_step > 0
          _num_steps = -1
        end
        return new(_order, _num_steps, _ds_step, radiation_damping_on, radiation_fluctuations_on)
      end
    end
  end
end

@def_integrator_struct(Yoshida) # Automatically selects split
@def_integrator_struct(MatrixKick)
@def_integrator_struct(BendKick)
@def_integrator_struct(SolenoidKick)
@def_integrator_struct(DriftKick)


# ========== Exact ===========================
struct Exact end

# ========== Explicit RK4 Tracking ==========
struct RungeKutta
  ds_step::Float64
  n_steps::Int
end

DEFAULT_RK4_DS_STEP = 0.2
function RungeKutta(; ds_step::Union{Float64, Nothing}=nothing, n_steps::Union{Int, Nothing}=nothing)
  # Get actual values (use provided or sentinel)
  _ds_step = ds_step === nothing ? -1.0 : ds_step
  _n_steps = n_steps === nothing ? -1 : n_steps
  
  # Error if both are explicitly set to positive values
  if _ds_step > 0 && _n_steps > 0
    error("Only one of ds_step or n_steps should be specified")
  end
  
  # If user sets n_steps (and it's positive), set ds_step to negative
  if _n_steps > 0
    return RungeKutta(-1.0, _n_steps)
  end
  
  # If user sets ds_step (and it's positive), set n_steps to -1
  if _ds_step > 0
    return RungeKutta(_ds_step, -1)
  end
  
  # Fallback: use defaults if both are negative/not set
  return RungeKutta(DEFAULT_RK4_DS_STEP, -1)
end