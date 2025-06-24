using Test,
      BeamTracking,
      Beamlines,
      JET,
      BenchmarkTools,
      GTPSA,
      StaticArrays

using BeamTracking: BunchView, KernelCall
BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
BenchmarkTools.DEFAULT_PARAMETERS.evals = 2

const D1 = Descriptor(6, 1)   # 6 variables 1st order
const D10 = Descriptor(6, 10) # 6 variables 10th order

function test_matrix(
  M_expected,
  kernel_call;
  type_stable=VERSION >= v"1.11", 
  no_scalar_allocs=!(any(t->eltype(t) <: TPS, kernel_call.args)), # only for non-parametric 
  rtol=nothing, 
  atol=nothing
)
  # Initialize bunch without spin
  v = transpose(@vars(D1))
  state = similar(v, State.T, 1)
  state .= State.Alive
  b = BunchView(state, v, nothing)

  # Set up kernel chain and launch!
  BeamTracking.launch!(b, kernel_call)

  # Set up tolerance kwargs
  kwargs = ()
  if !isnothing(atol)
    kwargs = pairs((;kwargs..., atol=atol))
  end
  if !isnothing(rtol)
    kwargs = pairs((;kwargs..., rtol=rtol))
  end

  # 1) Correctness
  @test isapprox(GTPSA.jacobian(b.v)[1:6,1:6], scalar.(M_expected); kwargs...)
  # 2) Type stability
  if type_stable
    @test_opt kernel_call.kernel(1, b, kernel_call.args...)
  end
  # 3) No scalar allocations
  if no_scalar_allocs
    v = [0.1 0.2 0.3 0.4 0.5 0.6]
    @test @ballocated(BeamTracking.launch!(b, $kernel_call; use_KA=false), 
    setup=(b = BunchView(copy($state), copy($v), nothing))) == 0
  end
end

function read_map(bmad_map_file::AbstractString)
  # Load reference data from file in isolated module to avoid polluting global namespace
  mod = Module()
  Base.include(mod, bmad_map_file)
  d_z = getfield(mod, :d_z)
  d_z == D10 || error("Please use a 10th order map for test_map")
  v_z = getfield(mod, :v_z)
  return v_z
end

function test_map(
  bmad_map_file::AbstractString,
  kernel_call;
  type_stable=VERSION >= v"1.11", 
  no_scalar_allocs=!(any(t->eltype(t) <: TPS, kernel_call.args)), # only for non-parametric 
  tol=1e-8
)
  v_expected = read_map(bmad_map_file)

  # Initialize bunch without spin
  v = transpose(@vars(D10))
  state = similar(v, State.T, 1)
  state .= State.Alive
  b = BunchView(state, v, nothing)

  # Set up kernel chain and launch!
  BeamTracking.launch!(b, kernel_call)

  # 1) Correctness
  @test coeffs_approx_equal(v_expected, b.v, tol)
  # 2) Type stability
  if type_stable
    @test_opt kernel_call.kernel(1, b, kernel_call.args...)
  end
  # 3) No scalar allocations
  if no_scalar_allocs
    v = [0.1 0.2 0.3 0.4 0.5 0.6]
    @test @ballocated(BeamTracking.launch!(b, $kernel_call; use_KA=false), 
    setup=(b = BunchView(copy($state), copy($v), nothing))) == 0
  end


  #= LineElement tracking test
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

  =#
end

#  Coefficient-wise approximate equality
function coeffs_approx_equal(v_expected, v_calculated, ϵ)
  n = GTPSA.numcoefs(v_expected[1])
  all_ok = true
  for i in 1:length(v_expected)
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
