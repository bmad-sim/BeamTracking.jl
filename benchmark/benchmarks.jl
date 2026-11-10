# Benchmarks of track!. Absolute timings depend on the machine, so compare revisions on the same
# machine, e.g. with AirspeedVelocity.jl (see .github/workflows/Benchmark.yml):
#
#   benchpkg BeamTracking --rev=main,dirty --script=benchmark/benchmarks.jl --add=Beamlines
#
# or locally with `include("benchmark/benchmarks.jl"); run(SUITE)`. The hardware-independent performance
# checks (allocations, type stability, timing ratios) are in test/performance_test.jl.
using BeamTracking, Beamlines, BenchmarkTools, Random

const SUITE = BenchmarkGroup()

const e = Species("electron")
# esr.jl defines the global `ring` itself, so the constant has another name
const esr = include(joinpath(@__DIR__, "..", "test", "lattices", "esr.jl"))
const ramp_ring = Beamline([deepcopy_no_beamline(ele) for ele in esr.line];
                           E_ref=Beamlines.R_to_E(e, esr.p_over_q_ref)*(1 + 1e3*Time()), species_ref=e)

# Elements with n_slices > 1 and the RK4 and Implicit tracking methods
quadratic_potential(x, y, s, t, k) = ((0.0, 0.0, 0.0, -k*x^2/2), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -k*x, 0.0, 0.0, 0.0))
const elements = Dict(
  "sliced quadrupole"  => LatElement(L=2.0, Kn1=0.1, tracking_method=MatrixKick(), n_steps=2, n_slices=4),
  "misaligned bend"    => LatElement(L=2.0, g=0.05, Kn0=0.05, tracking_method=BendKick(), n_steps=2, x_offset=1e-4, x_rot=1e-4, tilt=0.01),
  "rk4 quadrupole"     => LatElement(L=1.0, Kn1=0.2, tracking_method=RK4(), n_steps=10),
  "implicit potential" => LatElement(L=1.0, four_potential=quadratic_potential, four_potential_params=0.1,
                                     four_potential_normalized=true, tracking_method=Implicit(), n_steps=4),
)

# A Bunch with the reference species and energy of the Beamline (at t = 0) already set
function make_bunch(bl, N, T=Float64)
  Random.seed!(1)
  p_over_q_ref = bl.p_over_q_ref
  return Bunch(T.(1e-4 .* randn(N, 6)); species=e, p_over_q_ref=(p_over_q_ref isa Number ? p_over_q_ref : p_over_q_ref(0.0)))
end

# track! modifies the Bunch, so each sample tracks a new one
benchmark_track(bl, N, T=Float64; kw...) = @benchmarkable track!(b, $bl; $kw...) setup=(b = make_bunch($bl, $N, $T)) evals=1

SUITE["esr"] = BenchmarkGroup()
# 2, 10 and 15 particles: counts that do not fill the SIMD vectors
for N in (1, 2, 10, 15, 100, 1000)
  SUITE["esr"]["static", N] = benchmark_track(esr, N)
end
for N in (15, 1000)
  SUITE["esr"]["static Float32", N] = benchmark_track(esr, N, Float32)
end
for N in (1, 100)
  SUITE["esr"]["ramp", N] = benchmark_track(ramp_ring, N)
  SUITE["esr"]["ramp per particle", N] = benchmark_track(ramp_ring, N; ramp_update_each_particle=true)
end

# Prepared tracking (the kernel chains of the elements built once, as for the turns of a ring),
# only on revisions that have it (the script also runs on the base branch)
const BTBL = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)
if isdefined(BTBL, :prepare) && isdefined(BTBL, :PreKernelChain)
  SUITE["esr prepared"] = BenchmarkGroup()
  for T in (Float64, Float32), N in (1, 2, 10, 15, 100, 1000)
    chains = BTBL.prepare(make_bunch(esr, N, T), esr; notify=false)
    SUITE["esr prepared"][string(T), N] = @benchmarkable track!(b, $chains) setup=(b = make_bunch($esr, $N, $T)) evals=1
  end
end

SUITE["elements"] = BenchmarkGroup()
for (name, ele) in elements, N in (1, 1000)
  SUITE["elements"][name, N] = benchmark_track(Beamline([ele]; p_over_q_ref=-10.0, species_ref=e), N)
end
