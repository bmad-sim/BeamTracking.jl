# Runge-Kutta Tracking

`RungeKutta` tracks particles through static magnetic multipole fields with a
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

The RK body kernel supports drifts and static magnetic multipoles, including
solenoid, normal, and skew terms. It does not add support for RF or electric
fields, maps, patches, four-potentials, radiation, or spin tracking.

Fringe tracking is not implemented for RK. A bend with nonzero `e1` or `e2`
therefore produces an error instead of silently omitting the edge effect.

## Time-dependent values and reference ramping

`ramp_update_each_particle=true` uses the upstream per-particle reference-ramp
path. Time-dependent values are evaluated once for each particle, using that
particle's time at the element entrance. These evaluated values stay fixed
during every RK substep and during the `k1` through `k4` stages. They are not
reevaluated at intermediate RK positions.

## Callbacks

RK calls internal callbacks after each completed non-final substep. It does not
call them between the `k1`, `k2`, `k3`, and `k4` stages. The common tracking path
performs the final callback after element-exit processing.

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

The current low-level entry point is:

```julia
rk4_kernel!(i, coords, beta_0, tilde_m, charge, p0c, mc2,
            L, ds_step, n_steps, gx, gy, mm, kn, ks, p_over_q_ref)
```

`mm` contains the magnetic multipole orders. `kn` and `ks` contain the normal
and skew strengths. The kernel marks a particle as `STATE_LOST_PZ` when its
transverse velocity is unphysical.
