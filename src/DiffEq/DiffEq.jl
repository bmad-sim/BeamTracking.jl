module DiffEq
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
using ..DifferentialEquations
using ..DiffEqPhysics
export track!

# This module should be an extension for DifferantialEquations

Base.@kwdef struct Drift{T}
  L::T  # drift length [m]
end


function track!(bunch::Bunch, ele::DiffEq.Drift; work=get_work(bunch, Val{0}()))
  L = ele.L
  v = bunch.v
  beta_ref = sr_beta(bunch.beta_gamma_ref)
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)

  # Define Hamiltonian
  H(q, p, params) = p[3]/beta_ref - sqrt((p[3]+1/beta_ref)^2 - p[1]^2 - p[2]^2 - 1/bunch.beta_gamma_ref^2)

  prob = HamiltonianProblem(H, nothing, nothing, (0, L))
  function prob_func(prob, i, repeat)
    remake(prob, u0 = ArrayPartition([v.x[i], v.y[i], v.z[i]], [v.px[i], v.py[i], v.pz[i]]))
  end

  ensprob = EnsembleProblem(prob, prob_func=prob_func)
  sol = solve(ensprob, Yoshida6(), dt=L, trajectories=length(v.x), dense=false, save_end=true)
  for i in 1:length(v.x)
    v.x .= sol.u[i].u[end][1]
    v.y .= sol.u[i].u[end][2]
    v.z .= sol.u[i].u[end][3]
    v.px .= sol.u[i].u[end][4]
    v.py .= sol.u[i].u[end][5]
    v.pz .= sol.u[i].u[end][6]
  end

  return bunch
end 






end