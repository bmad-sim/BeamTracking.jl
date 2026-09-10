# Runge-Kutta Tracking

`RungeKutta` tracks particles through electromagnetic field sources with a
classical fourth-order Runge-Kutta (RK4) integrator. It uses the common tracking
path for reference-momentum updates, alignment, aperture checks, and callbacks.

## Configuration

Set either `ds_step` or `n_steps`:

- `ds_step` is the maximum requested step length in meters. The number of steps
  is `ceil(Int, L / ds_step)`, and the actual equal step length is `L / n_steps`.
- `n_steps` is the exact number of equal steps. The step length is
  `L / n_steps`.

Only one option can be positive. If neither option is set, the default is
`ds_step=0.2`. RK tracking requires a positive element length.

```julia
ele.tracking_method = RungeKutta()
ele.tracking_method = RungeKutta(ds_step=0.1)
ele.tracking_method = RungeKutta(n_steps=50)
```

## Field sources

A field source is a concrete callable object with the interface:

```julia
source(x, y, z, s) -> EMField
```

`EMField.E` is an `SVector` in V/m and `EMField.B` is an `SVector` in tesla.
RK tracking provides `ZeroField`, `MultipoleField`, `FunctionalField`, and
`SumField` source types.

- `ZeroField()` returns zero electric and magnetic fields.
- `MultipoleField(orders, normal, skew)` stores static multipole orders and
  non-integrated physical magnetic coefficients. Its field values are in
  tesla.
- `FunctionalField(evaluator, parameters)` stores the evaluator type and the
  parameter type in the source type. The evaluator receives
  `(x, y, z, s, parameters)` and returns an `EMField`.
- `FunctionalField(evaluator)` calls the evaluator with `(x, y, z, s)`.
- `SumField(sources...)` stores a tuple of concrete sources and evaluates the
  sum with static dispatch.

`field` sets the complete body field:

```julia
function uniform_field(x, y, z, s, parameters)
  v = zero(x)
  return EMField(
    v, v, v,
    v + parameters.Bx, v + parameters.By, v + parameters.Bz,
  )
end

source = FunctionalField(uniform_field, (Bx=0.0, By=0.1, Bz=0.0))
ele.tracking_method = RungeKutta(field=source, n_steps=20)
```

`additional_field` adds a source to the magnetic multipoles stored on the
element:

```julia
ele.tracking_method = RungeKutta(additional_field=source, n_steps=20)
```

The configured sources and their parameter types remain concrete in the RK
kernel.

### Parameter unpacking

Tuples and named tuples provide the standard parameter containers for a
`FunctionalField`. `SVector`s provide the coefficient containers for a
`MultipoleField`. Beamlines prepares these values before the tracking function
barrier:

1. `DefExpr` leaves are evaluated with the beamline `Context`.
2. `scalar_params=true` applies Beamlines scalarization to each parameter
   leaf.
3. `BatchParam` values are lowered and selected for each particle.
4. `TimeDependentParam` values are evaluated at each particle's element-entry
   time.

The field evaluator receives the resulting concrete parameter value. Ordinary
arrays can hold field-map data inside a functional source. Field source types
participate in `Adapt`, so backend array adaptation reaches nested source
parameters. A GPU evaluator uses GPU-compatible Julia operations and
device-compatible parameter storage.

### Custom sources

Custom callable objects can be passed directly to `field` or
`additional_field` when their fields are already concrete and do not require
parameter preparation:

```julia
struct UniformMagneticField{T}
  By::T
end

function (source::UniformMagneticField)(x, y, z, s)
  v = zero(x)
  return EMField(v, v, v, v, v + source.By, v)
end

ele.tracking_method = RungeKutta(field=UniformMagneticField(0.1))
```

Use `FunctionalField` when a custom source needs `DefExpr`, `BatchParam`,
`TimeDependentParam`, or scalarization. Keep these values in its tuple or named
tuple parameter container, following the same evaluator-and-parameters pattern
as Beamlines `MapParams`:

```julia
function uniform_magnetic_field(x, y, z, s, parameters)
  v = zero(x)
  return EMField(v, v, v, v, v + parameters.By, v)
end

source = FunctionalField(
  uniform_magnetic_field,
  (By=DefExpr{Float64}(c -> c.strength),),
)
ele.tracking_method = RungeKutta(field=source)
# The beamline Context must define strength.
```

Field preparation is explicit for `MultipoleField`, `FunctionalField`, and
`SumField`. It does not inspect arbitrary structs, array storage, or closure
captures. `FunctionalField` always preserves its evaluator and processes only
its parameters. Custom parameter structs are opaque; use tuples or named tuples
for values that require preparation. Custom types containing device storage
also need their own `Adapt` support.

For explicit SIMD, `x`, `y`, `z`, and `s` can be `SIMD.Vec` values. Field
evaluators can use `zero(x)` as a numeric carrier when combining coordinates
with scalar parameters, as in `uniform_field` above. The same evaluator also
works with floating-point, dual, and TPSA coordinates.

## Beamlines usage

The element must get its reference data from a `Beamline`, in the same way as
the upstream tracking methods.

```julia
using BeamTracking, Beamlines

species = Species("electron")
ele = Quadrupole(L=1.0, Kn1=0.1,
                 tracking_method=RungeKutta(ds_step=0.1))
line = Beamline([ele], p_over_q_ref=-3.0, species_ref=species)
bunch = Bunch(zeros(100, 6), p_over_q_ref=line.p_over_q_ref,
              species=species)

track!(bunch, line)
```

The RK body kernel tracks static electric and magnetic fields. Beamlines
magnetic multipoles are represented by `MultipoleField`, including solenoid,
normal, and skew terms. Bend body tracking uses reference curvature with zero
edge angles.

## Time-dependent values and reference ramping

`ramp_update_each_particle=true` uses the upstream per-particle reference-ramp
path. Time-dependent values are evaluated once for each particle, using that
particle's time at the element entrance. Each evaluated value stays fixed
during every RK substep and during the `k1` through `k4` stages.

## Callbacks

RK calls internal callbacks after each completed non-final substep. Each RK
step completes its `k1`, `k2`, `k3`, and `k4` stages before the callback. The
common tracking path performs the final callback after element-exit processing.

## Reference curvature

For a bend, reference curvature is split using `tilt_ref`:

```math
g_x = g_{ref}\cos(tilt_{ref}), \qquad
g_y = g_{ref}\sin(tilt_{ref}).
```

The path-length correction and geometric momentum terms are

```math
dh = g_x x + g_y y,
```

```math
\frac{dp_x}{ds} = \frac{F_x}{p_0}\frac{dt}{ds} + g_x\frac{p_s}{p_0},
\qquad
\frac{dp_y}{ds} = \frac{F_y}{p_0}\frac{dt}{ds} + g_y\frac{p_s}{p_0}.
```

## Low-level kernel

The low-level entry point is:

```julia
rk4_kernel!(i, coords, beta_0, tilde_m, charge, p0c, mc2,
            L, ds_step, n_steps, gx, gy, source)
```

`source` is a concrete callable field source. The kernel marks a particle as
`STATE_LOST_PZ` when its transverse velocity is unphysical.
