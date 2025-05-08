using Test,
      BeamTracking,
      Beamlines,
      JET,
      BenchmarkTools,
      GTPSA

BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
BenchmarkTools.DEFAULT_PARAMETERS.evals = 2

const D = Descriptor(6, 1)

function test_matrix(kernel, M_expected, args...; type_stable=VERSION >= v"1.11", no_allocs=true, tol=1e-14)
  n_temps = BeamTracking.MAX_TEMPS(parentmodule(kernel).TRACKING_METHOD())
  v = transpose(@vars(D))
  work = zeros(eltype(v), 1, n_temps)

  BeamTracking.launch!(kernel, v, work, args...)

  # 1) Correctness
  @test norm(GTPSA.jacobian(v)[1:6,1:6] - scalar.(M_expected)) < tol 
  if norm(GTPSA.jacobian(v)[1:6,1:6] - scalar.(M_expected)) > tol
    println("Actual: $(GTPSA.jacobian(v)[1:6,1:6])")
  end
  # 2) Type stability
  if type_stable
    @test_opt BeamTracking.launch!(kernel, v, work, args...)
  end
  # 3) No Allocations
  if no_allocs
    @test @ballocated(BeamTracking.launch!($kernel, $v, $work, $args...)) == 0 
  end
end

include("LinearTracking.jl")
include("BeamlinesExt.jl")
