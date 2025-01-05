module DiffEq
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

using ..GTPSA: @FastGTPSA!, GTPSA
using ..SciMLBase
using ..StaticArrays
using ..RecursiveArrayTools
using ..NaNMath


Base.@kwdef struct Drift{T, U}
  L::T  # drift length [m]
  solver::U
  ds::Float64 = L 
end

function track!(bunch::Bunch, ele::DiffEq.Drift) #; work=get_work(bunch, Val{0}()))
  L = ele.L
  solver = ele.solver
  ds = ele.ds
  particles = bunch.particles

  #= 
  Normalize to beta_gamma_ref = 0 for simplicity and speed
  We rewrite the fraction 

  1/[beta*(1 + p[3]^2 + 1/(beta*gamma)^2)] -> multiply top and bottom by gamma and move inside
  = gamma/(beta^2*gamma^2 + beta^2*gamma^2*p[3]^2 + 1) -> now evaluate right-sided limit gamma->1, beta->0
  = 1

  So for the dz/ds ODE this fraction goes to 1

  =#
  # Define Hamilton's equations:
  # SA = SArray = static array for no allocation
  pdot!(dp, p, q, params, t) = dp .= SA[0.0, 0.0, 0.0] # No change to momentum in drift
  qdot!(dq, p, q, params, t) = dq .= SA[p[1] / sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2), 
                                        p[2] / sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2),
                                        p[3] - (p[3] + 1)/sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2)]

  prob = DynamicalODEProblem(pdot!, qdot!, view(particles[1].v, 1:3), view(particles[1].v, 4:6), (0, L))

  function prob_func(prob, i, repeat)
    remake(prob, u0 = ArrayPartition(view(particles[i].v, 1:3), view(particles[i].v, 4:6)))
  end
  ensprob = EnsembleProblem(prob, prob_func=prob_func)
  sol = solve(ensprob, solver, dt = ds, trajectories=length(particles), dense=false, save_everystep=false, save_end=true)

  return sol
end 






end