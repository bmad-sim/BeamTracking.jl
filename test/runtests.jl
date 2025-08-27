using Test,
      BeamTracking,
      Beamlines,
      JET,
      BenchmarkTools,
      GTPSA,
      StaticArrays,
      ReferenceFrameRotations,
      SIMD

using BeamTracking: Coords, KernelCall, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST
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
  state = similar(v, UInt8, 1)
  state .= STATE_ALIVE
  coords = Coords(state, v, nothing)

  # Set up kernel chain and launch!
  BeamTracking.launch!(coords, kernel_call)

  # Set up tolerance kwargs
  kwargs = ()
  if !isnothing(atol)
    kwargs = pairs((;kwargs..., atol=atol))
  end
  if !isnothing(rtol)
    kwargs = pairs((;kwargs..., rtol=rtol))
  end

  # 1) Correctness
  @test isapprox(GTPSA.jacobian(coords.v)[1:6,1:6], scalar.(M_expected); kwargs...)
  # 2) Type stability
  if type_stable
    @test_opt kernel_call.kernel(1, coords, kernel_call.args...)
  end
  # 3) No scalar allocations
  if no_scalar_allocs
    v = repeat([0.1 0.2 0.3 0.4 0.5 0.6], 2)
    q = repeat([1.0 0.0 0.0 0.0], 2)
    state = [STATE_ALIVE STATE_ALIVE]
    @test @ballocated(BeamTracking.launch!(coords, $kernel_call; use_KA=false), 
    setup=(coords = Coords(copy($state), copy($v), copy($q)))) == 0
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

function read_spin_orbit_map(bmad_map_file::AbstractString)
  # Load reference data from file in isolated module to avoid polluting global namespace
  mod = Module()
  Base.include(mod, bmad_map_file)
  d_z = getfield(mod, :d_z)
  d_z == D10 || error("Please use a 10th order map for test_map")
  v_z = getfield(mod, :v_z)
  q_z = getfield(mod, :q_z)
  return (v_z, q_z)
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
  q = TPS64{D10}[1 0 0 0]
  state = similar(v, UInt8, 1)
  state .= STATE_ALIVE
  coords = Coords(state, v, q)

  # Set up kernel chain and launch!
  BeamTracking.launch!(coords, kernel_call)

  # 1) Correctness
  @test coeffs_approx_equal(v_expected, coords.v, tol)
  # 2) Type stability
  if type_stable
    @test_opt kernel_call.kernel(1, coords, kernel_call.args...)
  end
  # 3) No scalar allocations
  if no_scalar_allocs
    v = repeat([0.1 0.2 0.3 0.4 0.5 0.6], 2)
    q = repeat([1.0 0.0 0.0 0.0], 2)
    state = [STATE_ALIVE STATE_ALIVE]
    @test @ballocated(BeamTracking.launch!(coords, $kernel_call; use_KA=false), 
    setup=(coords = Coords(copy($state), copy($v), copy($q)))) == 0
  end


  #= LineElement tracking test
    if haskey(kwargs, :R_ref) && haskey(kwargs, :species)
      R_ref = kwargs[:R_ref]
    elseif haskey(kwargs, :E) && haskey(kwargs, :species)
      R_ref = BeamTracking.E_to_R(kwargs[:species], kwargs[:E])
    elseif haskey(kwargs, :p0c) && haskey(kwargs, :species)
      R_ref = BeamTracking.E_to_R(kwargs[:species], sqrt(kwargs[:p0c]^2 + BeamTracking.BeamTracking.massof(kwargs[:species])^2))
    else
      error("`R_ref`, `E` or `p0c`, as well as `species` must both be provided as keyword arguments")
    end
    
    if !haskey(kwargs, :ele)
      error("ele must be provided as a keyword argument")
    else
      coords = Bunch(v, species=kwargs[:species], R_ref=R_ref)
      v = track!(coords, kwargs[:ele]).v
      @test coeffs_approx_equal(v_expected, v, tol)
    end

  =#
end

# Coefficient-wise approximate equality
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


function quaternion_coeffs_approx_equal(q_expected, q_calculated, ϵ)
  sgn = ifelse(q_expected.q0[[0,0,0,0,0,0]] * q_calculated.q0[[0,0,0,0,0,0]] >= 0, 1, -1)
  components = (:q0, :q1, :q2, :q3)
  n = GTPSA.numcoefs(q_expected.q0)
  all_ok = true
    for cname in components
      v_expected = getfield(q_expected, cname)
      v_calculated = sgn * getfield(q_calculated, cname)
      for j in 0:n-1
          c1, c2 = v_expected[j], v_calculated[j]
          if abs(c1 - c2) > max(ϵ, ϵ * (abs(c1) + abs(c2)))
              println("Coefficients not equal: expected $cname[$j] = $c1, got $cname[$j] = $c2")
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
include("IntegrationTracking.jl")
include("BeamlinesExt.jl")