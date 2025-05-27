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
  # 2) Type stability
  if type_stable
    @test_opt BeamTracking.launch!(kernel, v, work, args...)
  end
  # 3) No Allocations
  if no_allocs
    @test @ballocated(BeamTracking.launch!($kernel, $v, $work, $args...)) == 0 
  end
end

function test_map(kernel, bmad_map_file::AbstractString, args...; type_stable=VERSION >= v"1.11", no_allocs=true, tol=1e-8)
  # Load reference data from file in isolated module to avoid polluting global namespace
  mod = Module()
  Base.include(mod, bmad_map_file)
  D = getfield(mod, :d_z)
  v_expected = getfield(mod, :v_z)

  # Construct test input variables
  n_temps = BeamTracking.MAX_TEMPS(parentmodule(kernel).TRACKING_METHOD())
  v = transpose(@vars(D))
  work = zeros(eltype(v), 1, n_temps)
  args = to_type.(args, TPS64{D})

  # Run kernel
  v = BeamTracking.launch!(kernel, v, work, args...)

  # 1) Correctness
  @test coeffs_approx_equal(v_expected, v, tol)
  # 2) Type stability
  if type_stable
      @test_opt BeamTracking.launch!(kernel, v, work, args...)
  end
  # 3) No allocations
  if no_allocs
      @test @ballocated(BeamTracking.launch!($kernel, $v, $work, $args...)) == 0
  end
end

#  Coefficient-wise approximate equality
function coeffs_approx_equal(v_expected, v_calculated, ϵ)
  n = GTPSA.numcoefs(v_expected[1])
  all_ok = true
  for i in 1:6
      for j in 0:n-1
          c1, c2 = v_expected[i][j], v_calculated[i][j]
          if abs(c1 - c2) > max(ϵ, ϵ * (abs(c1) + abs(c2)))
              println("Coefficients not equal: v_expected[$i][$j] = $c1, v_calculated[$i][$j] = $c2")
              println("Difference: $(abs(c1 - c2))")
              println("Tolerance:  $(max(ϵ, ϵ * (abs(c1) + abs(c2))))")
              all_ok = false
              break
          end
      end
      if !all_ok
          break
      end
  end
  return all_ok
end

function to_type(x, T)
  if isa(x, AbstractArray)
      return map(el -> convert(T, el), x)
  else
      return convert(T, x)
  end
end


include("LinearTracking.jl")
include("ExactTracking.jl")
include("BeamlinesExt.jl")
