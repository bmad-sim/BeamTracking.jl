# Working on BeamTracking.jl

BeamTracking-specific guide for contributors and AI coding agents: the architecture of the 
Beamlines extension, the physics conventions, the performance rules, how to verify changes and the 
current status. **The rules for how to work (mindset, environment, verification) are in 
`../SciBmad/AGENTS.md`**, which `CLAUDE.md` imports together with this file; they apply here 
strictly. Read both before changing code, and keep this file up to date.

BeamTracking.jl contains the tracking kernels and, in `ext/BeamTrackingBeamlinesExt/`, the code that 
tracks a `Bunch` through a [Beamlines.jl](https://github.com/bmad-sim/Beamlines.jl) `Beamline` or 
`LatElement`. Both packages are part of [SciBmad](https://github.com/bmad-sim/SciBmad.jl). Most work 
touches both repositories; their versions are coupled (BeamTracking 0.9 needs Beamlines 0.12).

## Environment

- Scratch environment: `Pkg.develop` `~/.julia/dev/BeamTracking`, `~/.julia/dev/Beamlines` (and 
  AtomicAndPhysicalConstants if checked out) and add the test dependencies (BenchmarkTools, JET, 
  GTPSA, StaticArrays, ReferenceFrameRotations, SIMD, KernelAbstractions, ForwardDiff, Distributions, 
  LinearAlgebra, Profile, Test). Run tests from it:
  `BEAMTRACKING_TEST_GROUP=<core|symplectic|batch_time> julia --project=<env> -e 'using Pkg; Pkg.test("BeamTracking")'`
  (no group: all of them). With unregistered Beamlines versions, `Pkg.test()` from the repo's own 
  project cannot resolve Beamlines.
- Single test files: run this script as `julia --project=<env> runsome.jl BeamlinesExt/beamlines_utils_test.jl …` 
  (paths relative to `test/`; it keeps the header of `test/runtests.jl`). A failing top-level 
  `@testset` stops the rest of its file, so for a negative control run the testset of interest alone.
  ```julia
  cd(joinpath(ENV["HOME"], ".julia/dev/BeamTracking/test"))
  src = read("runtests.jl", String)
  src = replace(src, "foreach(include, TEST_FILES)" => "foreach(f->include(joinpath(@__DIR__, f)), $(repr(ARGS)))")
  include_string(Main, src, joinpath(pwd(), "runtests.jl"))
  ```
- GPU (Apple silicon): an environment with Metal.jl; track a Float32 `Bunch(MtlArray(v))` and compare 
  with the same Bunch on the CPU; both should be about equally far from a Float64 run (Float32 
  rounding, e.g. ~1e-8 for one pass, ~1e-6 in z after 10 turns with RF). Test with `spin=true` too. 
  Metal has no Float64, so any Float64 constant reaching a kernel fails ("unsupported use of double 
  value"). Integer division is one: `1/2` and `sgn/2` are Float64. Write `x/2`, or for a standalone 
  constant a `@generated` function returning the literal in the coordinate precision, like 
  `c_light(T)` and `one_half(T)` in `src/utils/float32.jl` (`1//2` fails with SIMD `Vec`s). 
  To find the source, print `(e[1], e[2])` for `e in err.errors` of the `InvalidIRError` (the 
  message's own stack trace stops at "multiple call sites").
- CI (`.github/workflows/CI.yml`) runs Julia `lts` (1.10) and `1` (latest release) on Linux, macOS 
  and Windows, x64 and aarch64, sharded into the three test groups.
- Error message style: `"<Method> is not compatible with elements containing <XParams>. Try <Other> instead"`.

## Architecture of `ext/BeamTrackingBeamlinesExt/`

### Files

- `BeamTrackingBeamlinesExt.jl`: imports, the two `track!` methods (`track!(bunch, ::Beamline)` and 
  `track!(bunch, ::LatElement)`), includes.
- `utils_bl.jl`: `TrackReference`, `resolve_reference!`, `get_strengths`/`get_integrated_strengths`/
  `get_e_strengths`, `rf_omega_calc`, `rf_phi0_calc`, `rf_step_calc`, `drift_params(species, p_over_q_ref) -> (tilde_m, gamsqr_0, beta_0)`.
- `unpack_bl.jl`: `_SliceUniversalParams`, `LoweredElement`, `_track!` 
  (dynamic layer), `universal!` (after the function barrier).
- `prepared_bl.jl`: kernel chains (see "Kernel chains and prepared tracking"): `PreParameter`, 
  `PreKernelChain`, `process_preparameter(s)`, `run_chain!`, `run!`, `push_chain!`, `prepare`, 
  `track!(bunch, ::Vector{PreKernelChain})`. Not exported.
- `general_bl.jl`: `ibs_kick`, `alignment(kc, alignmentparams, bendparams, L, s, entering)`, `aperture`. 
  `ibs_kick` computes the bunch statistics and the number of particles from the live particles 
  only (`weights` zero for lost ones; `w`/`w_inv` there are the bend tilt rotations).
- `symplectic_bl.jl`: `fringe_in`/`fringe_out`, `integration_launcher(ker, params, photon_params, edge_params, tm, universalparams)` 
  (for `AbstractSymplectic` and `Exact`; picks `order_{two,…,ten}_integrator!` statically).
- One file per tracking method, each defining `fringe_body_fringe(tm, kc, bunch, ele)`: 
  `driftkick_bl.jl`, `matrixkick_bl.jl`, `solenoidkick_bl.jl`, `bendkick_bl.jl`, `implicit_bl.jl`, 
  `exact_bl.jl`, `rk4_bl.jl`. `sagancavity_bl.jl` has its own `universal!(…, tm::SaganCavity, …)` 
  with the same signature.

### Pipeline

1. `track!(bunch, bl)` resolves the reference once per call (`resolve_reference!` → 
   `TrackReference(species_ref, p_over_q_ref, field_scale)`), builds the `Val`s of the call-wide 
   flags once (`ramp_particle_energy_without_rf`, `ramp_update_each_particle`, `shift_time = bunch.species != species_ref`), 
   and calls `@noinline _track!` per element (a static call).
2. `_track!` is the **dynamic layer**:
   - Fetches all parameter groups in one pass over the element's `pdict` and then its parents' 
     (through `InheritParams`), child first, comparing types with `===`. This replaces ~15 
     `ele.X` lookups per element. A `ProtectParams` errors as in Beamlines' `getproperty`.
   - Applies `IgnoreParams` (`(ignoreparams::IgnoreParams).ignore_params`: **the assertion is 
     essential**, see "Performance rules"), `deval`s each group, replaces inactive groups and 
     `RFParams` with `rf_on = false` by `nothing`, optionally scalarizes (also `ref`).
   - `tm = universalparams.tracking_method` (error if `nothing`).
   - Calls `@noinline universal!(bunch, tm, ele, _SliceUniversalParams(universalparams, k), <groups 
     positionally>, ref, <Vals>, batch_start, sink)` once per slice (k = 1 of 1 for a whole element, 
     see "Slices").
3. `universal!` (after the barrier, all types known):
   - Evaluates `p_over_q_ref` at `bunch.t_ref` (a number; time dependent references are only kept 
     as functions for the per-particle ramp) and computes the reference times and energies at full 
     precision. **Only the `RefState` gets values lowered** to the coordinate number type 
     (`num_lower`); `push` lowers kernel arguments itself, and host-side differences (times, 
     energies, `bunch.t_ref`) must not be rounded to e.g. Float32.
   - Builds `KernelChain(Val{11}(), refstate)` and pushes, in order: IBS, entrance ramp 
     (`reference_momentum_shift!`), entrance aperture/alignment, `fringe_body_fringe`, exit 
     aperture/alignment, exit ramp (per-particle ramp only), `reference_time_shift!` (different 
     species only). The `LoweredElement` is built just before `fringe_body_fringe`.
   - Tracking (`sink` is the Bunch): `@noinline run_chain!(bunch, kc)` (evaluates the `PreParameter`s, 
     launches), then sets `bunch.t_ref` and `bunch.p_over_q_ref`. Otherwise `push_chain!(sink, …)` 
     (see "Kernel chains and prepared tracking").
   - `ref` (`TrackReference`) is the only source of the tracking `p_over_q_ref`. In 
     `fringe_body_fringe`, `ele.ref.p_over_q_ref` is that value evaluated at the element entrance, 
     except with the per-particle ramp, where it is the `TimeDependentParam` (so everything derived 
     from it, e.g. `drift_params` and `omega`, is evaluated per particle). SaganCavity gets the 
     unevaluated `ref`.

### `fringe_body_fringe` conventions

- Compatibility checks first (`isnothing(ele.xparams) || error(...)`), then `drift_params`.
- Each case ends in its own `return push(kc, integration_launcher(...))`: never merge different 
  `params` types into one call site.
- Pure drifts, bends and solenoids use the exact order-2 solution (`@set tm.order = 2`), **except 
  with radiation**: radiation is applied in the integrator steps, so the requested order sets how 
  accurately it is integrated. The condition constant-folds.
- Straight `edge_params`: `(a, tilde_m, Ksol, Kn0, Ks0, Kn1, Ks1, Kn2, Ks2)` (absent entries 
  `nothing`). Bend `edge_params`: `(a, tilde_m, beta_0, Kn0, w, w_inv, e1, e2)`. Edge angles are 
  part of the fringe, so they follow `fringe_at`.
- The straight hard-edge fringe (`fringe!` in `src/kernels/fringe.jl`, PTC's multipole fringe) 
  takes the field coefficients of `(x + i*y)^n`, i.e. `Kn/n!`: `Kn2/2` for sextupoles, in the 
  orbit and the spin (`bz` is the jump of the scalar potential, `(Kn/n!)/(n+1) Im[(x+i*y)^(n+1)]`). 
  It used `Kn2` like Bmad (whose `hard_multipole_edge_kick` uses `k2` for sextupoles but `Kn/n!` 
  for `sad_mult`): twice the fringe. Checked against Xsuite and MAD-NG (orbit, to 1e-18) and the 
  zero-edge-length limit of a Maxwell-consistent soft-edge sextupole with Thomas-BMT (orbit, 
  dominant spin rotation; tested in `BeamlinesExt_symplectic_test.jl`).
- Strengths: `get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)` with the 
  **element** length `L` (integrated strengths are per element length, also in a slice).
- RF: `omega = rf_omega_calc(ele.rfparams, ele.beamlineparams, ele.ref)`, `t_ref = (rf_phi0_calc(ele.rfparams, ele.ref.species_ref) - 3pi/2)/omega` (see "RF phase"), 
  cavity params `(radiation_params, beta_0, gamsqr_0, tilde_m, a, omega, t_ref, slip, E0_normalized, Ksol_or_nothing, Val{sol}(), mm, kn, ks)`. 
  Write the `Ksol`/`Val` pair as an `if`/`else` with literal `Val{true}()`/`Val{false}()`.
- `radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing`; `photon_params` 
  only if `isprimitivetype(eltype(bunch.coords.v)) && rad_fluctuations`.
- A `TimeDependentParam` compares equal to a number only if it is constant (`_isconst`), and any 
  operation with a non-constant one gives a non-constant one: once `Kn0` depends on time, so does 
  the rotated skew dipole `ks[1]`, even when zero. Checks for zero (e.g. the skew dipole in a bend, 
  `no_skew_dipole`) must then look at the stored values.
- `SolenoidKick` on an element without a solenoid prepends a zero solenoid component 
  (`pushfirst(mm, 0)`, `pushfirst(kn, zero(kn[1]))`).
- Implicit and the Exact electric dipole are wrapped in `bmad_to_mad!`/`mad_to_bmad!` with the 
  scalar potential taken as zero at the element ends; at interior slice ends they convert with 
  the potential (`callback_implicit!`/`callback_electric!`), as callbacks do.

### Kernels (BeamTracking `src/`)

- `KernelChain` holds a tuple of `KernelCall`s and a `RefState` (`t_enter, beta_gamma_enter, t_exit, beta_gamma_exit, L, g, ds_step, s0`). 
  `push` returns a new chain type. `make_kernel_call` time- and batch-lowers arguments; `_push` 
  `num_lower`s them to the coordinate type.
- The `order_*_integrator!`s only build the tuple of substep lengths of their scheme; 
  `integrate_steps!` (fringes, steps, radiation, callbacks) calls `ker(i, coords, s, params..., w)` 
  for each substep `w`, with `s` from the launch entrance. 
  Callbacks get `s0 + s` (`construct_main_callback`), i.e. from the element entrance.
- `compute_dt_ref` uses `find_tilde_m(ker)`, the position of the argument named `tilde_m`.
- Field kernels that depend on `s` take the element-relative offset: `rk4_kernel!(…, s0=0)` and 
  `implicit_integrator!(…, s0, L)`; their field time (`t = (s/beta_0 - z…)/c`) is computed from 
  that `s`, so it is element-local.
- `compute_time`: `t = -z/v + t_ref`. `cavity_kick!` takes the phase relative to the reference 
  particle at the kick (`t = t_ref - τ/c` after `bmad_to_mad!`).
- `reference_momentum_shift!`/`reference_time_shift!` use `muladd` (never a local `fma` that 
  shadows `Base.fma`).
- `sincu`/`sinhcu` (`src/utils/math_simd.jl`): Taylor series through `x^8` for `|x| < 0.1`, 
  `sin(x)/x` above, for every float type; within 1 eps (Float64, Float32) and ForwardDiff 
  derivatives to ~1e-14 (those of `sin(x)/x` cancel near 0). TPS arguments use GTPSA's own 
  `sincu`. A series to `x^2` below 0.004 had been off by up to 2e-12 (bends and coordinate 
  transforms with small angles; it made exact bends depend on `n_steps` at 1e-12).
- `TimeDependentParam{F<:TimeFunction}` is parametric, so operations on it build closures of known 
  types and evaluation is type stable (see "Performance rules").
- Float32/Float16 tracking: `num_lower(T, ::TimeFunction)` rebuilds the closure tree with the 
  captured Float64 numbers (lattice values, species constants) converted to `T` (`lower_captured`; 
  a closure's type parameters are its captured values' types, plus static parameters such as 
  `literal_pow`'s `N`), so time functions are evaluated in `T` (Metal has no Float64). It is never 
  called for Float64, and converting at evaluation instead (in the closures) would still put 
  Float64 constants into GPU kernels. `num_lower` probes the output type at `zero(T)`.

## Physics conventions

### Reference energy and species

- A Beamline has one `p_over_q_ref` (possibly a `TimeDependentParam`), the **signed** rigidity 
  R = pc/(c q), negative for electrons. It is stored in the first element's `InitialBeamlineParams`; 
  `trygetproperty(bl, :p_over_q_ref)` gives it or a `GetError`.
- Inheritance (`resolve_reference!`): an unset Bunch species/energy is taken from the Beamline, an 
  unset Beamline one from the Bunch; both unset is an error, and so is a Bunch `p_over_q_ref` 
  whose sign differs from its charge. If both energies are set but differ, the Bunch is ramped to 
  the Beamline reference at the first element.
- **A Bunch of a different species than `species_ref`** (e.g. positrons in an electron ring) sees 
  physically fixed fields, as in a real machine:
  - It is tracked at the same |Bρ| with the sign of its own charge: `TrackReference.p_over_q_ref` 
    = Beamline R × `field_scale` (±1). Normalized strengths are multiplied by `field_scale`, 
    unnormalized ones divided by the tracking R (`get_strengths`). For `field_scale == -1`, 
    normalized custom RK4 fields get the factor `-1.0` in place of `Val(true)` (`normalized_field_at`), 
    and normalized Implicit potentials go through the unnormalized path with `p_over_q_ref = -1.0` 
    (exact division; the kernels use `p_over_q_ref` for nothing else). No wrapper types.
  - Everything about the particle uses `bunch.species` (mass, charge, anomaly, radiation, IBS); 
    machine settings use `species_ref` (harmonic RF frequency, RF phase convention, design clock).
  - **Design clock:** `bunch.t_ref` is the time of the reference-species particle. The kernels 
    measure z against the Bunch reference particle; at the end of each element (slice) 
    `reference_time_shift!` re-references z to the design clock.
  - **RF slip:** in a cavity the Bunch reference drifts from the design particle by `s*slip`, 
    `slip = 1/v_bunch - 1/v_design` (0 for the same species); `cavity!` uses `t_ref + (s + L/2)*slip`.
  - SaganCavity: the Bunch reference changes by the same rigidity `dp_over_q_ref` as the design 
    particle (its own `dE_ref`); the kernels add the slip per drift from the design momenta passed 
    as `design_params = (mass_d, P0c_d, dE_d)` (`nothing` for the same species). Its `universal!` 
    stores the exit reference in `bunch.p_over_q_ref` after the launch, as the general one does 
    (it never did: with `dE_ref ≠ 0` the next Beamline then ramped the Bunch a second time).
  - `beta_0` cannot simply be changed instead: `tilde_m` must be the Bunch's, and the high precision 
    z formulas assume `beta_0 = 1/sqrt(1 + tilde_m^2)`.

### RF phase

- Every cavity method gives the energy gain `qV cos(rf_phi0_calc + omega*dt)`, `dt` the arrival time 
  after the reference particle (z < 0 is late): SaganCavity `cos(omega*(t_phi0 + dt))` with 
  `t_phi0 = rf_phi0_calc/omega`, `cavity_kick!` `sin(omega*(t_ref + dt))` with `t_ref` as above 
  (`-3pi/2` rather than `+pi/2` keeps `omega*t_ref` exactly 0 or -pi at the zero crossings with 
  `phi0 = 0`, as before). `rf_phi0_calc` adds `pi` for a negative `species_ref` and the 
  `zero_phase` offset: `Accelerating` 0, `AboveTransition` +pi/2, `BelowTransition` -pi/2.
- Physics check of the names (Beamlines' `PhaseRef` docstring): above transition a late particle 
  must lose energy. ESR (`AboveTransition`, above transition) has a stable one-turn longitudinal 
  map (trace 1.902) with both SolenoidKick and SaganCavity; with the zero crossings swapped it is 
  unstable (2.098). Tested in `beamlines_utils_test.jl` ("RF phase convention").
- Bmad's own elements differ in the time direction: lcavity `cos(2pi(phi0 + f dt))`, rfcavity 
  `sin(2pi(phi0 - f dt))`. `Accelerating` equals lcavity exactly; the zero crossings equal rfcavity 
  at `phi0 = 0` (the `bmad_maps/pure_rf.jl` check), with the opposite sign of `phi0` otherwise.
- Before `followups2`, SaganCavity had the zero crossings swapped and DriftKick/SolenoidKick 
  decelerated at `Accelerating`. The SaganCavity tests swap `AboveTransition`/`BelowTransition` 
  in their lattice instead, so they keep their numbers.

### Kernel chains and prepared tracking

- `universal!` builds the chain of an element (slice) and passes it to `sink`: the Bunch tracked 
  (`track!`: launched at once with `run_chain!`, then the Bunch's `t_ref` and `p_over_q_ref` are 
  set) or a `Vector` (`prepare`: stored as `PreKernelChain(kc, dt, p_over_q_ref_exit)`; the copy of 
  the Bunch used for building gets `p_over_q_ref_exit` unless it is only known at the call). 
  `track!(bunch, chains)` runs the stored chains: `run!` evaluates the `PreParameter`s, launches, 
  then `t_ref += dt` and sets `p_over_q_ref`. `universal!` is never part of a chain: every element 
  is lowered once, also with a time dependent reference energy and IBS. 
- `track!(bunch, bl)` deliberately does not go through `prepare` (the maintainer's decision): 
  storing the chains boxes each one and adds a dynamic call per element (+0.6 ms per ESR turn for 
  one particle). Code that tracks the same Beamline many times calls `prepare` itself (SciBmad's 
  `track` and `dynamic_aperture`; not `find_closed_orbit`, which tracks each coordinate type 
  about once).
- What depends on the Bunch at the call is a `PreParameter` in the chain (only when preparing; 
  `track!` computes everything at once): `bunch.t_ref` (`RefState` `t_enter`/`t_exit`, the 
  `reference_time_shift!` argument); with a time dependent reference energy the reference momentum 
  at the element, `PreParameter(b -> p(b.t_ref))`, and everything computed from it (the `RefState` 
  energies, `dt`, the exit reference, and without the per-particle ramp every strength computed in 
  the lowering); the Bunch's reference momentum (the ramp shifts); the IBS statistics 
  (`ibs_stats`, one `PreParameter`, so they are computed once per call: the IBS kernel takes them 
  as one tuple argument).
- `PreParameter` has arithmetic like `TimeDependentParam` (`+ - * / ^`, `literal_pow`, `sqrt`, 
  `sin`, …, `atan2`), so the lowering code is unchanged; comparisons error (a decision that 
  depends on the call cannot be made when preparing). `get_strengths`/`get_integrated_strengths` 
  keep a `PreParameter` reference momentum out of the strengths' type (they stay numbers; only 
  those divided by it are `PreParameter`s) and use the uniform branches as for 
  `TimeDependentParam`. `num_lower` of a `PreParameter` lowers its value at the call to the 
  coordinates' type; `process_preparameters` also looks inside tuples, named tuples and static 
  arrays of kernel arguments (only if the type contains a `PreParameter`: `has_preparameter` is a 
  `@generated` compile-time constant, else static prepared tracking was 12x slower).
- A `PreParameter` is one closure `bunch -> value` that calls the closures of the `PreParameter`s it 
  is computed from directly (like `TimeFunction`s). Evaluating nested `PreParameter`s through one 
  recursive method defeated type inference: a prepared ramp took 23.6 ms per ESR turn for one 
  particle (vs 2.2 ms rebuilding the chains) instead of 0.98 ms. Captured values must be concrete 
  (a captured type is stored as a `DataType`: runtime dispatch, caught by the JET test of 
  `universal!`).
- With a time dependent reference energy (without the per-particle ramp) the 
  `reference_momentum_shift!` is always in a prepared chain (whether the reference changes is only 
  known at the call); the kernel skips a zero shift (`(bg*pz - 0)/bg` rounded pz: 1 ulp in one of 
  the 38 bitwise cases, positrons through a SaganCavity, confirmed by removing the skip). 
  SaganCavity with a time dependent reference is not supported (also not unprepared) and errors 
  when preparing.
- Results (ESR, Julia 1.10, Apple M4, one particle, per turn): prepared static 0.50 ms, ramp 
  0.92-1.11 ms (2.17-2.21 rebuilding the chains at every call), per-particle ramp 1.63-1.70 ms 
  (5.87-6.10); 1000 particles unchanged (kernel bound). Prepared tracking is bitwise equal to 
  `track!` per turn for static lines, both ramp modes (Float64 and Float32), other species, 
  SaganCavity, IBS (also with a ramp), slices, spin. Negative controls: the reference energy 
  evaluated at preparation fails 10 checks, stale IBS statistics 1. Compile: first prepare + 
  prepared track! of a ramped ESR 6.3 s, as before.
- Decided against: the per-particle ramp's shift kernels reading the bunch reference and 
  `beta_gamma_exit` from the `RefState` (a new field); the `PreParameter` arguments already prepare 
  it fully. `ramp_update_each_particle` stays `false` by default (prepared per-particle ramp 1.7 vs 
  1.0 ms for one particle, 2.3x at 1000: the time function is evaluated per particle in the 
  kernels; the transforms for callbacks are disabled with it since d7d2882).
- Benchmarks: `SUITE["esr prepared"]` has both ramp modes (1 and 100 particles); 
  `test/performance_test.jl` checks prepared ramps allocate independently of the particle count 
  and take < 0.8x a static `track!` for one particle (0.25x and 0.42x on Apple M4).
- Reference momentum changes: the particles are shifted by exactly the change of the value the 
  Bunch stores (`convert(typeof(bunch.p_over_q_ref), …)`; shifted iff it changes; the per-particle 
  ramp ends at the stored value; SaganCavity's ramp the same). An approximate comparison 
  (`≈`, 1.5e-8 per element, 3.4e-4 for a Float32 reference) used to skip slow ramp steps while 
  storing them: a 3.3e-6 ramp over 1000 drifts was lost, every Float32 ramp below 3.4e-4 per 
  element was lost, and the Float32 per-particle ramp drifted (26% there). Test "Slow reference 
  momentum ramp" (`time_test.jl`; 6 of 12 fail without the fix); a Float32 ramp now agrees with 
  Float64 to 1.6e-8 (8.3e-6 before); everything else is bitwise unchanged (38 cases).

### Slices (`n_slices > 1`)

- Each slice is tracked like an element of length `L/n_slices` with its own `launch!`: the reference 
  is re-evaluated (ramped) and IBS recomputed at its start; the entrance fringe and aperture only 
  in the first slice, the exit ones only in the last; misalignments transform in at the slice 
  entrance and out at its exit, so the Bunch is in the lab frame between slices.
- `_SliceUniversalParams{A,F,RD,RF,FIRST,LAST}` (also for a whole element, k = 1 of 1; it replaced 
  `_ConcreteUniversalParams`) has the element `L`, per-slice `n_steps`/`ds_step`, `s0`, `L_slice`, 
  and `fringe_at` reduced to the element ends. Its `getproperty` returns `s0 = 0` (an `Int` 
  constant) for a first slice and `L_slice = L` for a whole element, as `_ConcreteUniversalParams` 
  did, so results and types are unchanged (38 tracking cases, ESR and sliced-FODO `twiss`/`track` 
  bitwise). Write the conditions as nested ternaries on the type parameters 
  (`key === :s0 ? (FIRST ? 0 : getfield(up, :s0)) : …`): with `key === :s0 && FIRST ? 0 : …` 
  `s0` was inferred `Union{Int, Float64}` and the `KernelChain` was boxed (+0.7 MB per ESR turn). 
  A whole element's struct is built by `_whole_element`, which branches on `fringe_at`, 
  `rad_damping` and `rad_fluctuations` (16 cases) so that each branch builds a type known at 
  compile time: built from those run-time values, the type is looked up in Julia's type cache at 
  every element, and the lookup grows with the number of type parameters (`_ConcreteUniversalParams` 
  124 ns per element, with `FIRST`/`LAST` 144-152 ns, `_whole_element` 24 ns). `track!` for one 
  particle on ESR: 3.78-3.82 ms with `_ConcreteUniversalParams`, 3.89-4.36 with the run-time 
  slice type, 3.16-3.20 ms now (allocations 3.36 -> 3.26 MB); compile time unchanged; results 
  bitwise equal (38 cases, ESR and sliced-FODO `twiss`/`track`).
- RF needs no `s0`: the kick phase is relative to the reference particle at the kick, and the 
  per-slice time shift re-aligns z with the design clock, so the slip restarts. Verified for all 
  `zero_phase`s, `phi0 ≠ 0`, DriftKick/SolenoidKick and a different species (sliced == pieces, 
  bitwise; the energy gain is independent of `n_slices`).
- Not supported: SaganCavity; elements with `PatchParams` or `MapParams`.
- The number of integrator steps (and of callbacks) of an element is always `n_slices*n_steps`: 
  elements tracked in one step (patches, maps, SaganCavity, whose kicks are set by `n_cells`) 
  error for `n_steps > 1`, and Beamlines forces `n_steps = n_slices = 1` for zero length. SciBmad's 
  `twiss` relies on this to know where each callback is.

### Misalignment geometry (verified)

- Straight elements: the body frame is shifted by `d = (x_offset, y_offset, z_offset)` and rotated 
  by `W = rot_quaternion(x_rot, y_rot, tilt)` = R_y(y_rot) R_x(x_rot) R_z(tilt) about the element 
  center. `track_alignment_straight_at_s!` works at any face `s`.
- Bends (Bmad convention): the rigid motion `x → C + d + W(x - C)` in the untilted frame at the 
  arc midpoint, `C` the chord midpoint (`coord_alignment_bend_mid`, `coord_alignment_bend_at_s`).
- Equivalent pieces: a half of a misaligned element (frame `P + Q y` at its own midpoint) has 
  rotation `Qᵀ W Q` and offsets `Qᵀ(C + d + W(P - C) - P) - C_h + Qᵀ W Q C_h` (straight: `Q = I`, 
  `C = 0`).

## Performance rules

- Transcendental functions (`sqrt`, `sin`, `cos`, `sinh`, `cosh`, …) of a SIMD `Vec` run lane by 
  lane (not vectorized) and are slow per particle too. Where the argument is small, a Taylor series 
  in Horner form (`x*(one(x)/n)` folds to a constant of the number type) behind `all(abs(x) < r)` 
  (one compare, plus one reduction for a `Vec`; a plain branch per thread on GPUs) is faster and as 
  accurate: `sincos_quaternion` (MatrixKick, spin) through x^10 for |x| < 1, 0.76 eps max vs 1.35 
  for the functions, derivatives 3e-16 vs 1e-8 (the old small-x branch). Gain: ESR prepared 
  1000 particles 217 -> 163 ms, the benchmark's 1-particle setting 1.80 -> 1.30 ms; compile time 
  unchanged; Float32 and Metal errors statistically unchanged (RMS over 100 particles). Test the 
  series against BigFloat with a negative control on a coefficient whose term is above round-off.
- The exact bend (`bend_kick.jl`) uses `asin_s`, `sincos_s`, `sin_s`, `atan2_s` (`math_simd.jl`), 
  series for floats and `Vec`s of them only (TPS and Duals keep the functions): ESR prepared 
  10-1000 particles 8-10% faster, one particle unchanged within noise (A/B with a series for `Vec`s 
  only: same); Float64 results change by 4e-15 after 2 ESR turns; Metal and CPU agree better 
  (same polynomial instead of two math libraries). Without 64-byte lanes the gain is about half 
  of the first round's (-15%), where the per-lane functions cost more.
- The series' branch on GPUs: a GPU group (warp) pays both paths only if its particles straddle 
  the threshold. The arguments are nearly the same for all particles: `k1 s^2/(4(1+δ))` in 
  MatrixKick (differs by δ ~ 1e-3), and in the exact bend px/pt (asin), the bend angle per step 
  plus ~px (sincos, spin) and ~x/L (atan2). A group splits only if `k1 s^2/4` is within ~0.1% of 
  1, or the bend angle per step within ~1e-3 of 0.1 rad. Measured on Metal (10^7 threads x 100 
  evaluations, Float32): `sincos_quaternion` uniform small arguments 2x faster, uniform large +4%, 
  alternating per thread (worst case, every group split) +20%; bend series with uniform arguments 
  10-20% faster than the functions (asin 13.8 vs 15.2 ms, sincos 14.8 vs 17.5, sin 12.9 vs 16.4, 
  atan2 13.7 vs 16.3). ESR on Metal is dominated by the ~4900 kernel launches per turn (10^4 and 
  10^5 particles take about as long): fusing launches would be a separate project.
- Apply the series to all floats, not only to SIMD `Vec`s: they are faster on every path. Without 
  explicit SIMD (`use_explicit_SIMD=false`, plain loop) ESR prepared 1000 particles 272-286 -> 
  207 ms (-25%), one particle 0.61 -> 0.53 ms: the compiler did not vectorize these kernels before 
  either (the calls to sin, asin, sqrt prevent it), and the series beats the scalar functions. With 
  SIMD 221-228 -> 149 ms (-33%, both series). SIMD lanes and the scalar loop do not give bitwise 
  equal results anyway (ESR 3 turns: 1.9e-13 Float64, 2.7e-5 Float32 without the series; the same 
  with them), so restricting the series to `Vec`s would not make them consistent either.
- Tried and rejected (no measurable gain, 2 A/B runs with and without SIMD): `sincu`/`sinhcu` 
  returning the series without evaluating `sin(x)/x` when all arguments are small 
  (`all(small) ? s : vifelse(...)`); the first round measured -3% before the bend series. TPS do 
  not use this `sincu`: BeamTracking extends `GTPSA.sincu`, whose `TPS` method is more specific.
- SIMD lane width (under evaluation, Nov 2026): 64-byte vectors (`max(native, 64 ÷ sizeof(T))`, 8 
  Float64 / 16 Float32) instead of the native width. Apple M4 (native 2 Float64, 4 Float32), ESR 
  prepared: Float64 1000 particles 146-151 -> 74 ms (2.0x), 100: 14.5-15.1 -> 8.2; Float32 1000: 
  67-68 -> 36 ms (1.86x). Alone it is slower when many particles are left for the scalar remainder 
  loop (Float32 4-15 and 31 particles up to 1.9x slower, Float64 2-7 up to 1.35x), so the remainder 
  first goes through native-width vectors: then every count is at least as fast as the native 
  width. Cost: a third compiled variant of every kernel chain (first `track!` on ESR 4.8-5.0 -> 
  6.1 s). x86 not measured: x86 Julia under Rosetta cannot load the arm64 GTPSA_jll override, 
  and Rosetta timings would not show native x86 anyway. Register budget: AVX2 holds as many 
  64-byte vectors (16 x 32 B) as NEON (32 x 16 B), AVX-512 is unchanged (native 64 B), SSE-only 
  CPUs hold half. Batching (`lane2vec`) works for any width.
  A/B script for another machine: `~/scibmad-benchmarks/lane-width/` (`run_lanes.jl`; with 
  `-C haswell` an AVX-512 machine runs AVX2 code, and HostCPUFeatures reports the AVX2 width 
  because it follows Julia's CPU target). GitHub's ubuntu runners (Benchmark.yml) are x86 too.
  Don't use `-C haswell` under WSL: it recompiles every package for that target and killed WSL 
  (and HostCPUFeatures warned it kept stale feature flags).
  Measured on x86 (the maintainer's laptop under WSL; Julia sees AVX2, native 4 Float64 / 8 
  Float32), best of 2 runs: Float64 1000 particles 151 -> 122 ms (1.24x), 100: 16.0 -> 13.7; 
  Float32 1000: 65 -> 48 ms (1.35x), 100: 9.8 -> 8.0; 1-31 particles equal within noise or faster 
  (largest slowdown +5% at Float32 9 particles; the two native runs differed by up to 20%).
  Further checks (Apple M4): results bitwise equal to the native width in all 38 bitwise cases 
  (vectors of any width compute lane by lane the same; only the scalar loop differs, which the 
  native-width remainder avoids); with 8 threads equal or faster at every count (10^4 particles 
  261-271 -> 158-161 ms); both full test suites pass. SciBmad precompile workload: 153 -> 158 s, 
  peak memory 10.7 -> 12.9 GB, binary 391 -> 459 MB (each kernel chain compiled for 3 variants). 
  Unrelated finding: `use_cpu_multithreading=true` costs ~42 ms per ESR turn below ~256 particles 
  (the `Threads.@threads` overhead per element; one particle single-threaded takes 0.5 ms).
  Portability: AVX-512 and SVE-512 are unchanged by construction (native = 64 bytes); not 
  measured: SVE-256 (Graviton 3, like AVX2), SSE-only x86 and the generic 16-byte targets (half the 
  register budget of NEON/AVX2 for these vectors).
- `benchmark/benchmarks.jl` has prepared tracking (`SUITE["esr prepared"]`, Float64 and Float32, 
  1/2/10/15/100/1000 particles; only on revisions that have `prepare`) and counts that do not 
  fill the SIMD vectors (2, 10, 15; Float32 15). `test/performance_test.jl` checks that prepared 
  tracking allocates about once per element independent of the particle count (4891 on ESR vs 
  81438 for `track!`) and that one particle takes < 0.4x `track!` (0.14x on Apple M4).
- `_track!` is dynamic; `universal!` is behind a `@noinline` function barrier. Anything fixed for 
  a whole `track!` call that changes types crosses the barrier as a `Val`, built once per call. 
  Never build `Val{runtime_bool}()`: use `cond ? Val{true}() : Val{false}()`.
- Per-element settings that change types cross the barrier as **type parameters of a struct that 
  already crosses it**, with a `getproperty` returning the type parameter, so call sites stay 
  unchanged and constant-fold: `fringe_at`, `rad_damping`, `rad_fluctuations` (`_SliceUniversalParams`), 
  `order`/`opt_scheme`/`use_newton` (the tracking methods in Beamlines), slice flags.
- In the dynamic layer every call on an untyped value is a dynamic dispatch, and its cost depends 
  on the method table (`in(::Type, ::Any)` was extremely slow). Type-assert values of known type, 
  and prefer `isnothing(x)` over calling functions on `nothing`.
- Never construct structs with many type parameters from untyped values, or use `@kwdef` keyword 
  constructors, in the dynamic layer: pass values positionally across the barrier.
- `ele.X` on a `LatElement` costs ~40 ns (it hashes `Type` keys, twice for `InheritParams` 
  children; every element in a `Beamline` is such a child). Don't use it per element in hot paths.
- A conditional `push` on a runtime value makes a union of `KernelChain` types; don't stack them. 
  Small `Union{Nothing,T}` are fine but multiply: 4–6 in one tuple is past union splitting.
- Both sides of ternaries/`ifelse` must have the same type (`zero(eltype(kn))`, not `0`). With 
  `TimeDependentParam`s they cannot (every operation is a new type): a broadcast `ifelse` over 
  per-order flags then builds its element type at runtime (~2 µs on 1.10, ~7 µs on 1.13 per 
  `get_strengths`), so `get_strengths` branches when all orders agree. In 
  kernels, values derived from the particle time are `SIMD.Vec`s: don't mix a scalar and a `Vec` 
  in one `ifelse` (a union that allocates per particle), and keep normalized strengths scalar.
- `TimeDependentParam{F}`: every operation gives a new concrete type. Use `<:TimeDependentParam` 
  in `::Type{…}` signatures; declare locals that get several such values as `TimeDependentParam` 
  (not `eltype(kn)`); build arrays with `T[T(x) for x in a]`, not `T.(a)` (broadcasting narrows the 
  element type).
- Kernels must not allocate, and must not loop over heterogeneous tuples with runtime indexing.
- **Compile time and code size.** The kernels (`@makekernel`) and the helpers they call are 
  `@inline`, so each is compiled into every call site, in every kernel chain type, for every 
  coordinate type (explicit SIMD and scalar paths, Float32, Dual, TPS with `@FastGTPSA`). Never 
  write the same kernel call several times in a row: loop over a tuple of what changes instead, so 
  there is one call site. The integrators used to write out every substep (up to 70 copies of the 
  kernel for order 10), which made the SciBmad precompile workload take 167 s and 12 GB of memory 
  and produced a 745 MB binary; with `integrate_steps!` (one call site, not inlined, so shared by 
  all chains) it is 72 s, 8.6 GB and 315 MB, the first `track!` through ESR compiles in 4.5 s 
  instead of 7.8 s (20 s → 4.6 s for a line with every tracking method), bitwise identical and 
  equally fast with explicit SIMD and on Metal.
- Removing `@inline` from the kernels themselves was tested (old PR #365 and again): the scalar 
  (no explicit SIMD) path becomes 2x slower and the explicit SIMD path ~7% slower, so the kernels 
  stay inlined. Without explicit SIMD the single call site costs ~10% on ESR (its alternating 
  quadrupoles; sextupoles and single quadrupoles got faster); writing the last step separately 
  only recovers part of it and costs compile time, so it was not done. `launch!` must stay 
  `@inline` (~1 ms per ESR turn otherwise).
- **Function arguments of non-inlined methods:** Julia does not specialize on a `Function` 
  argument that is only passed on or splatted (`ker(i, coords, params...)` is not a direct call); 
  it compiles the method for `::Function` and dispatches dynamically, here per particle (allocations 
  per particle, ~2x slower). Declare such arguments `ker::K` with `where {K}`. Inlined methods are 
  not affected, which is how it stays hidden until something is not inlined. The allocation test 
  in `test/performance_test.jl` (1 vs 100 particles, with and without explicit SIMD) catches it.
- **Julia 1.10 explicit SIMD codegen bug (x86-64 segfault, fixed in 1.11):** copying scalar fields 
  of an argument (e.g. the kernel chain, 8-aligned) into a tuple that contains a `Vec` gives the 
  load the destination's alignment (`emit_memcpy_llvm`, "memcpy_refined"); two merged `double` 
  copies become `movaps xmm, [arg+8]`, which faults when the argument is 16-aligned. aarch64 does 
  not fault on misaligned loads, so it only shows on x86 CI, and whether it faults depends on field 
  offsets (it can appear or vanish with unrelated changes). `_reg` (`src/time.jl`; the identity on 
  1.11+) rebuilds bits values from registers; the leaves of `teval` and `beval` go through it, so 
  per-particle time and batch evaluation are safe with explicit SIMD. **Any new code that builds a 
  tuple of per-particle `Vec`s together with scalars from kernel arguments must do the same.** 
  Base-only reproducer: `consume((x::NTuple{4,VecElement{Float64}}, args[2], args[3]))` in a 
  `@noinline` function with `args::NTuple{3,Float64}`; `code_native` on x86 shows 
  `movaps xmm0, [rdi + 8]`, and `movups` with the fields `bitcast` (as `_reg` does).
- KernelAbstractions 0.9 and 0.10: in 0.10 `GPU` is a deprecated alias of `Backend` and `CPU` is a 
  `Backend` (PoCL). Dispatch GPU paths on `::Backend` (the `::CPU` methods are more specific); a 
  "GPU" check is `!(backend isa CPU)`.
- `get_strengths` is `@inline` on purpose. Kernel setup runs once per element; it matters for many 
  elements with few particles (rings), kernels matter for many particles.

## Verification

### Correctness

- **Bitwise checks against the previous state:** a script that tracks a fixed bunch, e.g. 
  `v0 = [1e-3 1e-4 -1e-3 2e-4 1e-3 1e-4; -2e-3 0 1e-3 -1e-4 0 -2e-4]`, two turns through 
  ESR (`include("test/lattices/esr.jl")`) in Float64 and Float32 and as positrons; through a ramped 
  ring (`E_ref = E*(1 + a*Time())`) for both species with `ramp_update_each_particle` false and 
  true; and through short lines with time dependent strengths (uniform, mixed normalized/unnormalized, 
  mixed integrated, with DriftKick and MatrixKick); plus the cavities of 
  `test/lattices/sagan_cavity_lat.jl` and the batch tests when relevant. It `serialize`s every 
  `(coords.v, t_ref)` to a file given as argument. Run it before (`git stash`, or `git stash push -- <dir>`, 
  then `git stash pop`) and after, and compare with `isequal` per case. Commits by the maintainer 
  between rounds are the "before" state.
- **Equivalence tests** are the strongest physics checks, e.g. `test/BeamlinesExt/beamlines_slices_test.jl`: 
  an element with `n_slices = k` tracks bitwise like k elements of length `L/k`; a misaligned 
  element like adjusted halves; s/t-dependent fields like elements with shifted fields; a tilted 
  electric dipole with a rotated particle like an untilted one. Explain any non-bitwise residual 
  (e.g. differently rounded lever arms: a few ulps of `L`) and include a negative control.
- Tests in `test/` are grouped in `test/runtests.jl` (`TEST_GROUPS`); every `*_test.jl` file must 
  be in a group. Inside a `@testset`, nested functions share its locals, so don't reuse names.

### Performance

- Benchmark: ESR (`include("test/lattices/esr.jl")`, `ring`, 4,891 elements), 1 and 10 particles, 
  after a first `track!` so the Bunch has its reference. Timings drift ±0.5 ms between processes: 
  alternate A/B runs (`git stash`/file swaps) at least twice, use the minimum.
- Tools: split the cost into unpacking ("front") and chain building plus kernels ("back") by 
  redefining `universal!` with the same signature as a no-op; `Profile` flat output with `C=true`, 
  attributing `ijl_apply_generic` samples to the first Julia caller to find dynamic dispatch; JET 
  `@report_opt target_modules=(ext, BeamTracking) universal!(…)` with real arguments per element kind.
- `test/performance_test.jl` (core group) has the hardware-independent regression checks: 
  allocations of `track!` independent of the number of particles, allocation upper bounds per 
  Julia version (update them when allocations drop or for new Julia versions), JET `@test_opt` on 
  a drift's `universal!`, and timing ratios on the same machine. It cannot catch a uniform slowdown 
  with unchanged allocations; that needs a same-machine comparison of absolute timings:
- Compile cost (for changes to inlining, kernels or the kernel chain): the time of the first 
  `track!` in a fresh process (ESR, and a line with every tracking method), and the SciBmad 
  precompile workload with `precompile_workload = true` in a separate environment whose compiled 
  caches go to a scratch depot (`JULIA_DEPOT_PATH=<scratch depot>:~/.julia`, delete 
  `<scratch depot>/compiled` before each run), under `/usr/bin/time -l` (wall time, maximum 
  resident set size) plus the size of `compiled/v1.10/SciBmad/*.dylib`. Do not turn the workload 
  on in the environments used otherwise.
- `benchmark/benchmarks.jl` (`SUITE`, BenchmarkTools) and `.github/workflows/Benchmark.yml` 
  (AirspeedVelocity.jl: PR vs base on the same runner, results in the job summary). Locally: 
  `include("benchmark/benchmarks.jl"); run(SUITE)`.

## Beamlines facts

- Tracking methods live in Beamlines and hold only integrator settings: `DriftKick`, `MatrixKick`, 
  `BendKick`, `SolenoidKick` (`order`, `opt_scheme`), `Implicit` (also `use_newton`), `Exact`, 
  `SaganCavity`, `RK4`. There is no default: kind constructors choose one (Solenoid/RFCavity/
  CrabCavity → SolenoidKick, SBend → BendKick, Quadrupole → MatrixKick, Sextupole/Octupole/
  Multipole/Drift/kickers → DriftKick, Marker/BPM/Patch/Map and zero length → Exact, 
  BeamBeam → Implicit).
- Element settings in `UniversalParams`: `n_steps`/`ds_step` (per slice, one independent 
  variable), `n_slices`, `fringe_at::Fringe.T`, `rad_damping`, `rad_fluctuations`; IBS in 
  `IBSParams` (`ibs_model::IBSModel.T`). Setting `L = 0` switches to `Exact()`.
- Every element in a `Beamline` is an `InheritParams` child of the element it was made from. 
  Setting a property on a child writes to the shared parent; only a group put in the child's 
  `pdict` overrides.
- `BMultipoleParams{T,N}` stores `SizedVector`s with one eltype for all strengths and tilts, so 
  once one is a `DefExpr{T}`, the constant ones are wrapped as `DefExpr`s of the promoted type. 
  `DefExpr(f)` infers `T` with `Base.promote_op` (so `DefExpr(()->a)` with a non-`const` global `a`, 
  or `DefExpr(c->c.k)`, is a `DefExpr{Any}`, and makes them all `DefExpr{Any}`); with an explicit 
  `DefExpr{Float64}(…)` (or an inferable return type) they are all `DefExpr{Float64}`. `deval` of 
  a group with `DefExpr`s builds a new group (copying such arrays); of a group without any, it 
  returns the group itself (statically for concretely typed fields, by an `isa DefExpr` check for 
  untyped ones), so never mutate a `deval`ed group.
- `rawget(ele_or_group, key)` (exported, like `deval`) returns an independent variable as stored 
  (e.g. the `DefExpr`), defined per parameter group with a default for `AbstractParams`; dependent 
  ones (e.g. `Kn1L`/`Bn1` when `Kn1` was set last, `ds_step` when `n_steps` was set, `E_ref` when 
  `p_over_q_ref` was set, `s`) error and name the independent variable. A raw dependent quantity 
  would be stale as soon as another independent variable it uses changes (e.g. `p_over_q_ref` for `Bn1`).
- Energy conversions (`E_to_R`, `R_to_v`, `pc_to_R`, delta conversions, …) are in Beamlines 
  (unexported). BeamTracking only has `beta_gamma_to_v`.
- AtomicAndPhysicalConstants: `gyromagnetic_anomaly` is `NaN` for atomic species (spin tracking 
  gives NaN); `Species` equality compares every field.

## Status (Nov 2026)

The restructured extension (`n-slices`), Beamlines 0.12 and `followups2` (Float32 time dependence 
and batches, also on Metal; Dual/TPS element lengths; one RF phase convention for all cavity 
methods) are merged. Branch `noinline2`: integrators with one kernel call site (`integrate_steps!`, 
see "Performance rules"); it replaces the old `noinline` branch (PR #365).

ESR, Julia 1.10, Apple M4, per turn (Nov 2026, branch `prepare`): `track!` 1 particle 3.2 ms 
(3.26 MB allocated); prepared (`prepare`, as SciBmad's `track`/`dynamic_aperture` use it) 1/10/
100/1000 particles 0.50/1.3/8.0/75 ms (Float32 1000: 36 ms), ramp 1.0 ms, per-particle ramp 1.7 ms 
for one particle. Older: `main` 18.8 ms → 3.9 ms (`noinline2`)

Open items:
- SaganCavity does not support a ramped reference (`rf_step_calc` rounds a time dependent 
  `omega`, and the reference energy change between Beamlines would need ramp semantics).
- Further compile-size reductions (not done): every kernel chain is compiled twice for primitive 
  types (explicit SIMD loop and the scalar remainder loop); the remainder could use a masked SIMD 
  step. Of the 315 MB SciBmad binary, at most ~70 MB is tracking code (with no kernel inlined at all 
  it is 245 MB), the rest is SciBmad's own workload.
- Radiation fluctuations (`rad_fluctuations`) allocate 4 times per particle without explicit SIMD 
  (also on `main` before `noinline2`).
- Remaining unpacking costs (ESR profile, before `_whole_element`, which removed most of the first): the `_SliceUniversalParams` construction (dynamic 
  type with runtime type parameters) ~11%, data-dependent unions after the barrier ~6% (fixing 
  them must not slow down general multipoles), the barrier dispatch ~5%. `deval` of a parameter 
  group without any `DefExpr` returns the group itself (~4 ns; it used to rebuild it, 50 ns for 
  `UniversalParams` and array copies for `BMultipoleParams`).
- Per-particle ramp: mostly per-particle evaluation of time-function chains. Making all strengths 
  one common `TimeDependentParam` type halves setup but makes normalized strengths `Vec`s in the 
  kernels (30% slower for 1000 particles), so it was not done.
- Possibly: error in `resolve_reference!` for spin tracking of species with a `NaN` anomaly.
- The pre-restructuring code differs from the current one by ~2e-10 after one ESR turn from a 
  nonzero start: only in the RF cavities, which it integrated at order 2 (`SciBmadStandard`); 
  `RFCavity` now defaults to `SolenoidKick` (order 4). With `DriftKick(order=2)` the results are 
  bitwise equal; order 4 is converged to ~1e-17 per cavity, order 2 was off by ~1e-11.
