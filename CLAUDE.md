# Working on BeamTracking.jl

Guide for contributors and AI coding agents. Read all of it before changing code. It covers how to 
work here, the architecture of the Beamlines extension, the physics conventions, the performance 
rules, how to verify changes, and the current status. Keep it up to date: when a session settles 
something (a convention, a pitfall, a measured result, an open question), add it here.

BeamTracking.jl contains the tracking kernels and, in `ext/BeamTrackingBeamlinesExt/`, the code that 
tracks a `Bunch` through a [Beamlines.jl](https://github.com/bmad-sim/Beamlines.jl) `Beamline` or 
`LatElement`. Both packages are part of [SciBmad](https://github.com/bmad-sim/SciBmad.jl). Most work 
touches both repositories; their versions are coupled (BeamTracking 0.9 needs Beamlines 0.12).

## How to work here

### Mindset

- **Physical correctness and compiler-friendly code are the top priorities**, in that order, then 
  minimal diffs.
- **Minimal changes.** Keep the existing structure and the shape of functions you edit. Don't add 
  helper functions, abstractions or reorganize code unless it is necessary or asked for. Match the 
  surrounding style and comment density. 2-space indentation everywhere.
- **Verify, don't assume.** Re-derive physics from the actual kernel and Beamlines source, not from 
  memory. Never guess a field, property or function name: look it up (e.g. it is `n_steps`, not 
  `num_steps`). When asked to check something, say plainly what you verified, what you could not, 
  and what was wrong before.
- **Measure, don't guess** (performance). Every change is measured A/B against the previous state; 
  keep only changes with a clear gain and revert the rest. Report what each change bought.
- **Every change is checked against the previous results**, bitwise when it should not change 
  them (see "Verification"). When results change, explain exactly why (e.g. a fused multiply-add 
  rounds differently) and quantify it.
- **Find root causes.** A test that fails after a refactor may be an outdated test or a real bug; 
  find out which before changing either. If a change would affect physics results or behavior 
  beyond the task, describe the issue and the reasoning to the maintainer before making it.
- **Use controls.** An equivalence test is only meaningful if a plausible mistake breaks it: also 
  run a negative control (e.g. the wrong rotation convention, `s` restarting per slice) and show 
  that it fails by a large margin.
- **Don't commit or push**; the maintainer commits between rounds. Check `git log`/`git status` 
  before diffing so a new commit is not mistaken for lost work. Don't run the whole test suite 
  for every small change; run the relevant files and bitwise checks, and the full suite at the end.
- **Report concisely**: what changed, what it bought, what was verified and how, what is left. 
  Flag anything surprising found along the way (bugs elsewhere, performance oddities) without 
  necessarily fixing it.

### Environment

- Develop with a scratch environment that `Pkg.develop`s local clones of BeamTracking, Beamlines 
  and AtomicAndPhysicalConstants and adds the test dependencies (BenchmarkTools, JET, GTPSA, 
  StaticArrays, ReferenceFrameRotations, SIMD, KernelAbstractions, ForwardDiff, Distributions, 
  Test). Run tests from that environment:
  `BEAMTRACKING_TEST_GROUP=<core|symplectic|batch_time> julia --project=<env> -e 'using Pkg; Pkg.test("BeamTracking")'`.
  (With unregistered Beamlines versions, `Pkg.test()` from the repo's own project cannot resolve 
  Beamlines.) To run single files, include `test/runtests.jl`'s header with a different file list.
- CI (`.github/workflows/CI.yml`) runs Julia `lts` (1.10) and `1` (latest release) on Linux, macOS 
  and Windows, x64 and aarch64, sharded into the three test groups. Check changes on both Julia 
  versions when allocations, inference or timings matter (`juliaup`, e.g. `julia +release`).
- Reference source for the upstream packages: `git clone --depth 1` of `bmad-sim/BeamTracking.jl`, 
  `bmad-sim/Beamlines.jl`, `bmad-sim/AtomicAndPhysicalConstants.jl`. Local branches may be ahead.
- Error message style: `"<Method> is not compatible with elements containing <XParams>. Try <Other> instead"`.

## Architecture of `ext/BeamTrackingBeamlinesExt/`

### Files

- `BeamTrackingBeamlinesExt.jl`: imports, the two `track!` methods (`track!(bunch, ::Beamline)` and 
  `track!(bunch, ::LatElement)`), includes.
- `utils_bl.jl`: `TrackReference`, `resolve_reference!`, `get_strengths`/`get_integrated_strengths`/
  `get_e_strengths`, `rf_omega_calc`, `rf_phi0_calc`, `rf_step_calc`, `drift_params(species, p_over_q_ref) -> (tilde_m, gamsqr_0, beta_0)`.
- `unpack_bl.jl`: `_ConcreteUniversalParams`, `_SliceUniversalParams`, `LoweredElement`, `_track!` 
  (dynamic layer), `universal!` (after the function barrier).
- `general_bl.jl`: `ibs_kick`, `alignment(kc, alignmentparams, bendparams, L, s, entering)`, `aperture`.
- `symplectic_bl.jl`: `fringe_in`/`fringe_out`, `integration_launcher(ker, params, photon_params, edge_params, tm, universalparams)` 
  (for `AbstractSymplectic` and `Exact`; picks `order_{two,…,ten}_integrator!` statically).
- One file per tracking method, each defining `fringe_body_fringe(tm, kc, p_over_q_ref, bunch, ele)`: 
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
   - With `n_slices == 1`, calls `@noinline universal!(bunch, tm, ele, _ConcreteUniversalParams(universalparams), <groups positionally>, ref, <Vals>, batch_start)`; 
     otherwise calls it once per slice with `_SliceUniversalParams(universalparams, k)` (see "Slices").
3. `universal!` (after the barrier, all types known) builds `LoweredElement` from the groups, then:
   - Evaluates `p_over_q_ref` at `bunch.t_ref` (a number; time dependent references are only kept 
     as functions for the per-particle ramp) and computes the reference times and energies at full 
     precision. **Only the `RefState` gets values lowered** to the coordinate number type 
     (`num_lower`); `push` lowers kernel arguments itself, and host-side differences (times, 
     energies, `bunch.t_ref`) must not be rounded to e.g. Float32.
   - Builds `KernelChain(Val{11}(), refstate)` and pushes, in order: IBS, entrance ramp 
     (`reference_momentum_shift!`), entrance aperture/alignment, `fringe_body_fringe`, exit 
     aperture/alignment, exit ramp (per-particle ramp only), `reference_time_shift!` (different 
     species only).
   - `@noinline launch!(bunch.coords, kc)`, then sets `bunch.t_ref` and `bunch.p_over_q_ref`.
   - `ref` (`TrackReference`) is the only source of the tracking `p_over_q_ref`. The `p_over_q_ref` 
     argument of `fringe_body_fringe` is that value evaluated at the element entrance.

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
- Strengths: `get_strengths(ele.bmultipoleparams, L, p_over_q_ref, ele.ref.field_scale)` with the 
  **element** length `L` (integrated strengths are per element length, also in a slice).
- RF: `omega = rf_omega_calc(ele.rfparams, ele.beamlineparams, ele.ref)`, `t_ref = (rf_phi0_calc(ele.rfparams, ele.ref.species_ref) - pi/2)/omega`, 
  cavity params `(radiation_params, beta_0, gamsqr_0, tilde_m, a, omega, t_ref, slip, E0_normalized, Ksol_or_nothing, Val{sol}(), mm, kn, ks)`. 
  Write the `Ksol`/`Val` pair as an `if`/`else` with literal `Val{true}()`/`Val{false}()`.
- `radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing`; `photon_params` 
  only if `isprimitivetype(eltype(bunch.coords.v)) && rad_fluctuations`.
- `SolenoidKick` on an element without a solenoid prepends a zero solenoid component 
  (`pushfirst(mm, 0)`, `pushfirst(kn, zero(kn[1]))`).
- Implicit and the Exact electric dipole are wrapped in `bmad_to_mad!`/`mad_to_bmad!` with the 
  scalar potential taken as zero at the element ends; at interior slice ends they convert with 
  the potential (`callback_implicit!`/`callback_electric!`), as callbacks do.

### Kernels (BeamTracking `src/`)

- `KernelChain` holds a tuple of `KernelCall`s and a `RefState` (`t_enter, beta_gamma_enter, t_exit, beta_gamma_exit, L, g, ds_step, s0`). 
  `push` returns a new chain type. `make_kernel_call` time- and batch-lowers arguments; `_push` 
  `num_lower`s them to the coordinate type.
- Integrators call `ker(i, coords, s, params..., ds_step)` with `s` from the launch entrance. 
  Callbacks get `s0 + s` (`construct_main_callback`), i.e. from the element entrance.
- `compute_dt_ref` uses `find_tilde_m(ker)`, the position of the argument named `tilde_m`.
- Field kernels that depend on `s` take the element-relative offset: `rk4_kernel!(…, s0=0)` and 
  `implicit_integrator!(…, s0, L)`; their field time (`t = (s/beta_0 - z…)/c`) is computed from 
  that `s`, so it is element-local.
- `compute_time`: `t = -z/v + t_ref`. `cavity_kick!` takes the phase relative to the reference 
  particle at the kick (`t = t_ref - τ/c` after `bmad_to_mad!`).
- `reference_momentum_shift!`/`reference_time_shift!` use `muladd` (never a local `fma` that 
  shadows `Base.fma`).
- `TimeDependentParam{F<:TimeFunction}` is parametric, so operations on it build closures of known 
  types and evaluation is type stable (see "Performance rules").

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
    as `design_params = (mass_d, P0c_d, dE_d)` (`nothing` for the same species).
  - `beta_0` cannot simply be changed instead: `tilde_m` must be the Bunch's, and the high precision 
    z formulas assume `beta_0 = 1/sqrt(1 + tilde_m^2)`.

### Slices (`n_slices > 1`)

- Each slice is tracked like an element of length `L/n_slices` with its own `launch!`: the reference 
  is re-evaluated (ramped) and IBS recomputed at its start; the entrance fringe and aperture only 
  in the first slice, the exit ones only in the last; misalignments transform in at the slice 
  entrance and out at its exit, so the Bunch is in the lab frame between slices.
- `_SliceUniversalParams{A,F,RD,RF,FIRST,LAST}` has the element `L`, per-slice `n_steps`/`ds_step`, 
  `s0`, `L_slice`, and `fringe_at` reduced to the element ends. `_ConcreteUniversalParams`'s 
  `getproperty` returns the constants `s0 = 0`, `L_slice = L`, `first_slice = last_slice = true`, 
  so the `n_slices = 1` code is unchanged.
- RF needs no `s0`: the per-slice time shift re-aligns z with the design clock, so the slip restarts.
- Not supported: SaganCavity; elements with `PatchParams` or `MapParams`.

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

- `_track!` is dynamic; `universal!` is behind a `@noinline` function barrier. Anything fixed for 
  a whole `track!` call that changes types crosses the barrier as a `Val`, built once per call. 
  Never build `Val{runtime_bool}()`: use `cond ? Val{true}() : Val{false}()`.
- Per-element settings that change types cross the barrier as **type parameters of a struct that 
  already crosses it**, with a `getproperty` returning the type parameter, so call sites stay 
  unchanged and constant-fold: `fringe_at`, `rad_damping`, `rad_fluctuations` (`_ConcreteUniversalParams`), 
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

- **Bitwise checks against the previous state:** before a change, track a fixed bunch (e.g. 
  `Bunch([1e-3 1e-4 -1e-3 2e-4 1e-3 1e-4; -2e-3 0 1e-3 -1e-4 0 -2e-4])` two ESR turns, plus 
  Float32, a different species, a ramped ring) and `serialize` `coords.v` and `t_ref`; compare 
  bitwise after. Use `git stash`/copies of single files to get the "before" results.
- **Equivalence tests** are the strongest physics checks, e.g. `test/BeamlinesExt/beamlines_slices_test.jl`: 
  an element with `n_slices = k` tracks bitwise like k elements of length `L/k`; a misaligned 
  element like adjusted halves; s/t-dependent fields like elements with shifted fields; a tilted 
  electric dipole with a rotated particle like an untilted one. Explain any non-bitwise residual 
  (e.g. differently rounded lever arms: a few ulps of `L`) and include a negative control.
- Compare against `main` (a worktree plus its own environment) for behavior that should not change.
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
- `BMultipoleParams{T,N}` stores `SizedVector`s; `deval` copies them.
- Energy conversions (`E_to_R`, `R_to_v`, `pc_to_R`, delta conversions, …) are in Beamlines 
  (unexported). BeamTracking only has `beta_gamma_to_v`.
- AtomicAndPhysicalConstants: `gyromagnetic_anomaly` is `NaN` for atomic species (spin tracking 
  gives NaN); `Species` equality compares every field.

## Status (Nov 2026)

Open PRs: BeamTracking `n-slices` → `main` (restructured extension, species model, performance, 
`n_slices`, performance tests and benchmarks) and Beamlines `newuniversal` → `main` (Beamlines 
0.12). The BeamTracking PR needs Beamlines 0.12 registered.

ESR, Julia 1.10, Apple silicon, 1 particle: `main` 18.8 ms → 4.6 ms (103,961 allocations); ramped 
reference 4.7 ms; per-particle ramp 1/100/1000 particles 9/57/483 ms (Julia 1.13: 9.5 ms for 1).

Open items:
- The Julia 1.10 x86 SIMD segfault is fixed by `_reg` (see "Performance rules"). Two tests still 
  skip explicit SIMD on 1.10 x86 and can probably be enabled: `batch_simd_supported` in 
  `test/RungeKuttaTracking_test.jl` and `test/BeamlinesExt/beamlines_ibs_test.jl`.
- GPU/Float32 gaps (all also on `main`): `zero`/`one(::BatchParam)` use `first(b.batch)`, scalar 
  indexing that errors for GPU arrays (e.g. `MatrixKick` with a `BatchParam(MtlArray)` strength; 
  `DriftKick` works); a Float32 `TimeDependentParam` strength with `MatrixKick` and the Float32 
  per-particle ramp (even with a Float32 `p_over_q_ref` ramp) become Float64 time functions and 
  fail `num_lower`. Float32 batch gathers with explicit SIMD hit an upstream SIMD.jl bug.
- Remaining unpacking costs (ESR profile): `deval(::UniversalParams)` ~16%, the 
  `_ConcreteUniversalParams` construction (dynamic type with runtime type parameters) ~11%, 
  `deval(::BMultipoleParams)` array copies ~8%, data-dependent unions after the barrier ~6% (fixing 
  them must not slow down general multipoles), the barrier dispatch ~5%.
- Per-particle ramp: mostly per-particle evaluation of time-function chains. Making all strengths 
  one common `TimeDependentParam` type halves setup but makes normalized strengths `Vec`s in the 
  kernels (30% slower for 1000 particles), so it was not done.
- Possibly: error in `resolve_reference!` for spin tracking of species with a `NaN` anomaly; make 
  `rf_omega_calc` use the evaluated `p_over_q_ref` (behavior change, deferred); 
  `_ConcreteUniversalParams.ds_step::Float64` would fail for a TPS or Dual `L`.
- `main` and the restructured code differ by ~4e-11 in z after one ESR turn from a nonzero start; 
  predates the performance work, not investigated.
