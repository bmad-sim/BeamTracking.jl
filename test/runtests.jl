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

function test_map(kernel, bmad_map_file::AbstractString, args...; type_stable=VERSION >= v"1.11", kernel_test=true, TPS_params=false, no_allocs=true, tol=1e-8, kwargs...)
  # Load reference data from file in isolated module to avoid polluting global namespace
  mod = Module()
  Base.include(mod, bmad_map_file)
  D = getfield(mod, :d_z)
  v_expected = getfield(mod, :v_z)

  # Construct test input variables
  v = transpose(@vars(D))

  if kernel_test # Kernel tracking test
    n_temps = BeamTracking.MAX_TEMPS(parentmodule(kernel).TRACKING_METHOD())
    work = zeros(eltype(v), 1, n_temps)
    if TPS_params
      args = map(el -> el === nothing ? nothing : convert.(TPS64{D}, el), args)
    end

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

  else # LineElement tracking test
    if haskey(kwargs, :Brho_ref) && haskey(kwargs, :species)
      Brho_ref = kwargs[:Brho_ref]
    elseif haskey(kwargs, :E) && haskey(kwargs, :species)
      Brho_ref = BeamTracking.calc_Brho(kwargs[:species], kwargs[:E])
    elseif haskey(kwargs, :p0c) && haskey(kwargs, :species)
      Brho_ref = BeamTracking.calc_Brho(kwargs[:species], sqrt(kwargs[:p0c]^2 + BeamTracking.massof(kwargs[:species])^2))
    else
      error("`Brho_ref`, `E` or `p0c`, as well as `species` must both be provided as keyword arguments")
    end
    
    if !haskey(kwargs, :ele)
      error("ele must be provided as a keyword argument")
    else
      b = Bunch(v, species=kwargs[:species], Brho_ref=Brho_ref)
      v = track!(b, kwargs[:ele]).v
      @test coeffs_approx_equal(v_expected, v, tol)
    end
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


include("LinearTracking.jl")
include("ExactTracking.jl")
include("BeamlinesExt.jl")
