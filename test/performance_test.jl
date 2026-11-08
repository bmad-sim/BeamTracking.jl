# Performance regression tests that do not depend on the hardware: allocation counts (which only
# depend on the Julia version), type stability, and timing ratios measured on the same machine.
# Absolute timings are compared between revisions on the same machine with benchmark/benchmarks.jl.
using Random

@testset "Performance" begin
  ring = include("lattices/esr.jl")
  e = Species("electron")
  E0 = Beamlines.R_to_E(e, ring.p_over_q_ref)
  ramp_ring = Beamline([deepcopy_no_beamline(ele) for ele in ring.line]; E_ref=E0*(1 + 1e3*Time()), species_ref=e)
  sliced_line = Beamline([
    LatElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(), n_slices=3),
    LatElement(L=2.0, g=0.05, Kn0=0.05, tracking_method=BendKick(), n_slices=2, x_offset=1e-4, tilt=0.01),
    LatElement(L=1.0, Kn1=0.2, tracking_method=RK4(), n_steps=4, n_slices=2),
  ]; p_over_q_ref=-10.0, species_ref=e)
  function bunch(bl, N)
    Random.seed!(1)
    return Bunch(1e-4 .* randn(N, 6); species=e, p_over_q_ref=(bl === sliced_line ? -10.0 : ring.p_over_q_ref))
  end
  function allocations(bl, N; kw...)
    b = bunch(bl, N)
    track!(b, bl; notify=false, kw...)
    track!(b, bl; kw...)
    return @allocations track!(b, bl; kw...)
  end
  # Coverage and forced bounds checking change the generated code
  instrumented = Base.JLOptions().code_coverage != 0 || Base.JLOptions().check_bounds == 1

  # The kernels do not allocate, so the allocations of track! (setting up the kernels of each
  # element) do not depend on the number of particles
  cases = (("static", ring, (;)), ("static without explicit SIMD", ring, (; use_explicit_SIMD=false)), ("ramp", ramp_ring, (;)),
           ("ramp per particle", ramp_ring, (; ramp_update_each_particle=true)), ("slices", sliced_line, (;)))
  allocs = Dict(name => allocations(bl, 1; kw...) for (name, bl, kw) in cases)
  for (name, bl, kw) in cases
    @test allocations(bl, 100; kw...) == allocs[name]
  end

  # Upper bounds, ~10% above the allocations on Julia 1.10.13 and 1.13.3 for the ESR (4891 elements).
  # Lower them when the allocations are reduced, and update them for new Julia versions
  if !instrumented
    bounds = VERSION < v"1.11" ? Dict("static" => 115_000, "ramp" => 123_000, "ramp per particle" => 143_000) :
                                 Dict("static" => 124_000, "ramp" => 132_000, "ramp per particle" => 152_000)
    for (name, bound) in bounds
      @test allocs[name] <= bound
    end
  end

  # No dynamic dispatch after the function barrier for a drift
  ext = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
  bl = Beamline([Drift(L=1.0)]; p_over_q_ref=-10.0, species_ref=e)
  b = bunch(sliced_line, 4)
  track!(b, bl; notify=false)
  ele = bl.line[1]
  up = Beamlines.deval(ele.UniversalParams, bl.context)
  blp = Beamlines.deval(ele.BeamlineParams, bl.context)
  ref = ext.TrackReference(e, -10.0, 1)
  @test_opt target_modules=(ext, BeamTracking) ext.universal!(b, ele.tracking_method, ele, ext._ConcreteUniversalParams(up), nothing, nothing, nothing,
                                                               nothing, nothing, nothing, blp, nothing, nothing, nothing, nothing, nothing, ref, Val(false), Val(false), Val(false), 1)

  # Timing ratios on the same machine (with generous margins): ramping the reference energy costs
  # about the same as tracking with a static reference (~1.0 on Julia 1.10 and 1.13), and the per
  # particle ramp ~2x for one particle
  if !instrumented
    function min_time(bl; kw...)
      t = @benchmark(track!(bb, $bl; $kw...), setup=(bb = $bunch($bl, 1); track!(bb, $bl; notify=false, $kw...)), evals=1, samples=20, seconds=10)
      return minimum(t).time
    end
    t_static = min_time(ring)
    @test min_time(ramp_ring) / t_static < 1.5
    @test min_time(ramp_ring; ramp_update_each_particle=true) / t_static < 4
  end
end
