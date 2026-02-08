# ========== Fringe  ===========================
abstract type AbstractFringeAt end
struct NoEnd <: AbstractFringeAt end
struct BothEnds <: AbstractFringeAt end
struct EntranceEnd <: AbstractFringeAt end
struct ExitEnd <: AbstractFringeAt end


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
      fringe_at::AbstractFringeAt
      ibs_num_particles::Int
  
      function $(esc(name))(; order::Int=4, num_steps::Int=-1, ds_step::Float64=-1.0, radiation_damping_on::Bool=false, radiation_fluctuations_on::Bool=false, fringe_at::AbstractFringeAt=BothEnds(), ibs_num_particles::Int=0)
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
        return new(_order, _num_steps, _ds_step, radiation_damping_on, radiation_fluctuations_on, fringe_at, ibs_num_particles)
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
struct Exact
  fringe_at::AbstractFringeAt
  ibs_num_particles::Int
  
  function Exact(; fringe_at::AbstractFringeAt=BothEnds(), ibs_num_particles::Int=0)
    return new(fringe_at, ibs_num_particles)
  end
end