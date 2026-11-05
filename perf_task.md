# Task: recover (and beat) `track!` performance on the ESR lattice

Read `CLAUDE.md` first: it holds all the context (code layout, how I want you to work, the reference/species model, type stability rules).

## The regression

Single particle, ESR ring (`include("test/lattices/esr.jl")`, `ring`), Julia 1.10.13, Apple silicon:

| | Time | Allocations | Memory |
|---|---|---|---|
| BeamTracking + Beamlines `main` | 14.3 ms | 216,441 | 6.60 MiB |
| New structure (BeamTracking `v0.9.0-1` + Beamlines `newuniversal`) | 20.7 ms | 221,659 | 9.63 MiB |

10 particles: 16.5 ms on `main` vs 23.1 ms on the new structure.

Branches: the new structure is **BeamTracking `v0.9.0-1` together with Beamlines `newuniversal`**. The baseline is `main` on both. Always switch the two repos as a pair.

Benchmark (keep it identical on both branches; the second `track!` is so the bunch already has species and energy set):

```julia
using BeamTracking, Beamlines, BenchmarkTools
include("test/lattices/esr.jl")
b0 = Bunch(zeros(1,6)); track!(b0, ring)
@btime track!($b0, $ring)
b10 = Bunch(zeros(10,6)); track!(b10, ring)
@btime track!($b10, $ring)
```

## What the numbers say

- The ring has 4,891 elements. The new branch does about **+5,200 allocations (about 1 more per element)**, **+3 MiB (about 650 bytes more per element)**, and **+6.4 ms (about 1.3 µs more per element)**.
- Going from 1 to 10 particles costs about the same on both branches, so the regression is in the **per-element setup** (unpacking, building the kernel chain, dispatch), not in the kernels that run per particle.
- The physics agrees: after one turn `t_ref` is identical, and the coordinates match `main` to about 1e-24 absolute. Any change must keep that.
- `main`'s `track!(bunch, bl)` calls the keyword `track!(bunch, ele; ...)` for every element, which should be *slower* than the new direct `_track!` call. So the cost is elsewhere.

## Suspects (unverified, ranked)

1. **The kernel chain loses its concrete type after `fringe_body_fringe`.** `integration_launcher` picks `order` (5 values), `fin`/`fout` (`Val`s from the runtime `fringe_at`) and `opt_scheme` at run time, so its return type can be a union of dozens of `KernelCall` types. That is past Julia's union-splitting limit, so `kc` would be inferred as abstract. Every later `push` (exit aperture/alignment, ramp, time shift) and `launch!` would then be dynamic dispatches that box values. Check with `@code_warntype` or Cthulhu on `universal!` for a lowered Drift and a lowered Quadrupole. If this is it, pass `fringe_at`, `order` and `opt_scheme` across the `universal!` barrier as `Val`s, the same way as the ramp flags and `shift_time`. Also compare with how `main` built the drift/quadrupole kernel calls.
2. **The cost of the barrier call itself.** `universal!` is called dynamically for each element with three runtime `Val`s, keyword arguments, and a larger `LoweredElement`. It now has 14 fields, including `ref::TrackReference`, whose `Species` holds a `String`. Measure the extra allocation with `Profile.Allocs`. Things to try:
   - drop `kwargs...` on that path
   - keep only `species_ref` and `field_scale` in `LoweredElement`, not the whole `TrackReference`
   - compute the `Val`s once per `track!` call instead of once per element
3. **Unpacking added in the restructure, if `main` doesn't do it:** `deval` of every param group, `isactive`, `IgnoreParams` handling, `scalarize`, `_ConcreteUniversalParams`, `LoweredElement` construction.
4. **`KernelChain` grew from `Val{10}` to `Val{11}`.** This is probably minor, but cheap to rule out.

## Plan

1. **Reproduce** both numbers with the script above: `main`/`main` vs `v0.9.0-1`/`newuniversal`. Beamlines changed too (`LatElement`, tracking methods moved there), so switch both repos together.
2. **Bisect** `git log main..v0.9.0-1` in BeamTracking, together with `git log main..newuniversal` in Beamlines (matching commits by date), to find which change or changes cost the time and the allocations.
3. **Profile:**
   - `Profile.Allocs.@profile` on `track!(b0, ring)` to find where the extra allocation per element comes from
   - `@profile` or `@profview` for the time
   - JET `@report_opt` / `@code_warntype` on `_track!`, `universal!` and `fringe_body_fringe` with real lowered elements of each kind in the ring (Drift, Quadrupole, SBend, Sextupole, RF cavity, Marker)
4. **Experiment freely**, but measure every change. Keep the minimal-change style from `CLAUDE.md`, and don't add helper functions unless one clearly pays for itself.
5. **Done means:**
   - time and allocations at or below `main` for 1 and 10 particles
   - results unchanged against `main` (bit-for-bit, or a stated tolerance)
   - the full test suite passes

   Report what each change bought.