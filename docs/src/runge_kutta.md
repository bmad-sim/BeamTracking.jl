# Runge-Kutta Tracking

The `RungeKutta` tracking method provides a 4th-order Runge-Kutta (RK4) integration scheme for tracking particles through arbitrary electromagnetic fields. This method is particularly useful for elements with complex fields that cannot be handled analytically.

## Overview

The Runge-Kutta tracking implementation uses a classical 4th-order Runge-Kutta method to numerically integrate the relativistic equations of motion for charged particles in electromagnetic fields. The implementation is optimized for GPU/SIMD compatibility using branchless operations and stack-allocated StaticArrays.

## Configuration

The `RungeKutta` tracking method can be configured with the following parameters:

### Parameters

- **`ds_step`**: The step size for integration (in meters). If not specified, defaults to `0.2` meters when neither parameter is set.

- **`n_steps`** : The number of integration steps. If specified and positive, the step size is calculated as `h = L / n_steps` where `L` is the element length.

**Important**: Only one of `ds_step` or `n_steps` should be specified. If both are set to positive values, an error will be raised. If neither is specified, `ds_step=0.2` is used by default.

### Examples

```julia
# Use default step size (0.2 meters)
ele.tracking_method = RungeKutta()

# Specify step size explicitly
ele.tracking_method = RungeKutta(ds_step=0.1)

# Specify number of steps
ele.tracking_method = RungeKutta(n_steps=50)
```

## Usage

### Basic Usage with Beamlines.jl

```julia
using BeamTracking, Beamlines

# Create an element with Runge-Kutta tracking
ele = Quadrupole(L=1.0, k1=0.1)
ele.tracking_method = RungeKutta(ds_step=0.1)
bl = Beamline([ele], R_ref=1e6)

# Create a bunch and track
bunch = Bunch(zeros(100, 6), R_ref=1e6, species=Species("electron"))
track!(bunch, ele)
```
The Runge-Kutta method works with all thick elements that have BMultipoleParams. Other thick elements are treated as drifts.

### Low-Level Kernel Interface

For advanced use cases, the `rk4_kernel!` function provides direct access to the integration routine:

```julia
rk4_kernel!(i, coords, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2,
            s_span, n_steps, g_bend, mm, kn, ks)
```

**Parameters:**
- `i`: Particle index
- `coords`: Particle coordinate array
- `beta_0`, `gamsqr_0`, `tilde_m`: Reference particle parameters
- `charge`: Particle charge (in units of elementary charge)
- `p0c`, `mc2`: Reference momentum and rest mass energy
- `s_span`: Integration range `(s_start, s_end)`
- `n_steps`: Number of integration steps
- `g_bend`: Curvature parameter (1/ρ for bends, 0 otherwise)
- `mm`: StaticArray of multipole orders
- `kn`: Normal multipole strengths
- `ks`: Skew multipole strengths

The multipole arrays define the magnetic field configuration. For example, a quadrupole with `k1=0.5` would use `mm=SA[2]`, `kn=SA[0.5]`, `ks=SA[0.0]`.

## Physics

The Runge-Kutta method integrates the relativistic equations of motion expressed in terms of arc length `s` as the independent variable. The state vector is:

```math
\mathbf{u} = [x, p_x, y, p_y, z, p_z].
```

The derivatives `du/ds` are calculated from the relativistic electromagnetism:

### Position Derivatives

```math
\frac{dx}{ds} = v_x \frac{dt}{ds}, \quad \frac{dy}{ds} = v_y \frac{dt}{ds}
```

where `v_x, v_y` are the transverse velocity components and `dt/ds` accounts for the relationship between arc length and time.

### Momentum Derivatives

```math
\frac{dp_x}{ds} = \frac{F_x}{p₀c} \frac{dt}{ds} + g_{bend} \frac{p_z}{p₀}, \quad \frac{dp_y}{ds} = \frac{F_y}{p₀c} \frac{dt}{ds}
```

where `F_x, F_y` are the Lorentz force components and `g_bend` is the curvature parameter (nonzero only in bend elements).

### Longitudinal Coordinate

```math
\frac{dz}{ds} = \text{rel\_dir} \left(\frac{\beta}{\beta₀} - 1\right) + \text{rel\_dir} \frac{\sqrt{1-v_t²} - 1 - dh_{bend}}{\sqrt{1-v_t²}} + \frac{d\beta}{ds} \frac{z}{\beta}
```
At the moment `rel_dir` is hardcoded to be 1. In the future this could be extended to track particles going backwards.

### Energy Deviation

```math
\frac{dp_z}{ds} = \frac{dp}{ds} / p₀c
```

where `dp/ds` is calculated from the work done by the electromagnetic fields.

## Implementation Details

### RK4 Algorithm

The 4th-order Runge-Kutta method uses the standard four-stage scheme:

1. **k₁**: Evaluate derivatives at current state `u(s)`
2. **k₂**: Evaluate derivatives at `u(s) + (h/2) k₁`
3. **k₃**: Evaluate derivatives at `u(s) + (h/2) k₂`
4. **k₄**: Evaluate derivatives at `u(s) + h k₃`

The final update is:

```math
u(s+h) = u(s) + \frac{h}{6}(k₁ + 2k₂ + 2k₃ + k₄)
```


### Particle Loss Condition

The implementation includes detection of unphysical particle states:

- **Unphysical momenta**: If the transverse velocity squared `v_t² ≥ 1`, the particle is marked as lost (`STATE_LOST_PZ`) and won't be tracked.

### Multipole Field Calculation

Electromagnetic fields are computed using the `multipole_em_field` function, which calculates field components from magnetic multipole parameters:

```julia
multipole_em_field(x, y, z, s, mm, kn, ks) -> (Ex, Ey, Ez, Bx, By, Bz)
```

**Parameters:**
- **`mm`**: StaticArray of magnetic multipole orders (0=solenoid, 1=dipole, 2=quadrupole, 3=sextupole, etc.)
- **`kn`**: Normal multipole strengths (normalized and not integrated)
- **`ks`**: Skew multipole strengths (normalized and not integrated)

**Field computation:**
- For `m=0` (solenoid): Returns longitudinal field `Bz`
- For `m≥1` (dipole, quadrupole, etc.): Computes transverse fields `Bx`, `By` using a Horner-like scheme for efficient polynomial evaluation
