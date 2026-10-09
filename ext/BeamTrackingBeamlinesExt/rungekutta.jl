# RungeKutta only constructs the body field integration kernel, which includes the hard edge fringes.
# Unpacking, reference-ramp, alignment, aperture, and callback
# are handled by the shared unpacking step.

# Each supported field group contributes parallel tuples of evaluators,
# parameter payloads, and units flags. Generic kernel parameter preparation
# handles these tuples without field-specific lowering or adaptation methods.
@inline function runge_kutta_field(bmultipoleparams, L, p_over_q_ref)
  !isactive(bmultipoleparams) && return ((), (), ())
  mm = getfield(bmultipoleparams, :order)
  bn, bs = get_strengths(bmultipoleparams, L, p_over_q_ref)
  parameters = mm isa Integer ? (SA[mm], SA[bn], SA[bs]) : (mm, bn, bs)
  return ((BeamTracking.multipole_field,), (parameters,), (Val(true),))
end

# Hard edge fringe parameters for the fringe! kernels, built as for the symplectic methods.
# A bend gets the curved (Hwang) fringe using its dipole and quadrupole strengths. Otherwise
# the straight fringe is used for the solenoid, dipole, quadrupole, and (if
# tm.multipole_fringe_on) higher order multipole components. Custom field functions get no
# fringe kick. A nonzero bend edge angle is only allowed if the fringe at that end is on.
# Also returns the solenoid strength for the canonical/mechanical momentum shift at the edges.
@inline function runge_kutta_edge_params(tm::RungeKutta, bunch, bendparams, bmultipoleparams, L, p_over_q_ref)
  fin = fringe_in(tm.fringe_at)
  fout = fringe_out(tm.fringe_at)

  if isactive(bendparams)
    (bendparams.e1 == 0 || fin isa Val{true}) ||
      error("RungeKutta tracking with a nonzero bend edge angle e1 requires the entrance fringe to be on")
    (bendparams.e2 == 0 || fout isa Val{true}) ||
      error("RungeKutta tracking with a nonzero bend edge angle e2 requires the exit fringe to be on")
  end

  tm.fringe_at == Fringe.NoEnd && return nothing, nothing, fin, fout

  tilde_m, _, _ = BeamTracking.drift_params(bunch.species, p_over_q_ref)
  a = gyromagnetic_anomaly(bunch.species)

  if isactive(bmultipoleparams)
    mm = getfield(bmultipoleparams, :order)
    kn, ks = get_strengths(bmultipoleparams, L, p_over_q_ref)
    if mm isa Integer
      mm, kn, ks = SA[mm], SA[kn], SA[ks]
    end
  else
    mm, kn, ks = nothing, nothing, nothing
  end

  if isactive(bendparams)
    Kn0 = zero(L)
    Kn1 = nothing
    if !isnothing(mm)
      for j in 1:length(mm)
        if mm[j] == 1
          Kn0 = kn[j]
          ks[j] ≈ 0 || error("A skew dipole field cannot yet be used with a bend fringe")
        elseif mm[j] == 2
          Kn1 = kn[j]
        end
      end
    end
    ntilt = -bendparams.tilt_ref
    if ntilt ≈ 0
      w = nothing
      w_inv = nothing
    else
      w = rot_quaternion(0, 0, ntilt)
      w_inv = inv_rot_quaternion(0, 0, ntilt)
    end
    edge_params = (a, tilde_m, Kn0, Kn1, w, w_inv, bendparams.e1, bendparams.e2,
                   bendparams.edge1_int, bendparams.edge2_int)
    return edge_params, nothing, fin, fout
  end

  isnothing(mm) && return nothing, nothing, fin, fout

  Ksol = nothing
  Kn0 = nothing
  tilt0 = 0
  for j in 1:length(mm)
    if mm[j] == 0
      Ksol = kn[j]
    elseif mm[j] == 1
      Kn0 = sqrt(kn[j]^2 + ks[j]^2)
      tilt0 = atan2(ks[j], kn[j])
    end
  end
  if tilt0 ≈ 0
    w0 = nothing
    w0_inv = nothing
  else
    w0 = rot_quaternion(0, 0, tilt0)
    w0_inv = inv_rot_quaternion(0, 0, tilt0)
  end
  mm_f, kn_f, ks_f = fringe_multipoles(tm, mm, kn, ks)
  if isnothing(Ksol) && isnothing(Kn0) && isnothing(mm_f)
    return nothing, nothing, fin, fout
  end
  return (a, tilde_m, Ksol, Kn0, w0, w0_inv, mm_f, kn_f, ks_f), Ksol, fin, fout
end

@inline runge_kutta_custom_field(::Nothing) = ((), (), ())

@inline function runge_kutta_custom_field(params::EMFieldParams)
  isnothing(params.em_field) && return ((), (), ())
  return ((params.em_field,), (params.em_field_params,),
          (Val{params.em_field_normalized}(),))
end

@inline function runge_kutta_body(
  tm::RungeKutta,
  kc,
  p_over_q_ref,
  bunch,
  bendparams,
  bmultipoleparams,
  patchparams,
  rfparams,
  mapparams,
  fourpotentialparams,
  emultipoleparams,
  em_field_params,
  L,
)
  isnothing(bunch.coords.q) || error("RungeKutta tracking does not support spin tracking")
  L > 0 || error("RungeKutta tracking requires a positive element length")
  !isactive(patchparams) || error("RungeKutta tracking does not support patch elements")
  !isactive(rfparams) || error("RungeKutta tracking does not support RF fields")
  !isactive(mapparams) || error("RungeKutta tracking does not support map elements")
  !isactive(fourpotentialparams) || error("RungeKutta tracking does not support FourPotentialParams")
  !isactive(emultipoleparams) || error("RungeKutta tracking does not support electric multipoles")

  edge_params, edge_ksol, fin, fout = runge_kutta_edge_params(tm, bunch, bendparams, bmultipoleparams, L, p_over_q_ref)

  if isactive(bendparams)
    g_ref = bendparams.g_ref
    tilt_ref = bendparams.tilt_ref
    gx = g_ref * cos(tilt_ref)
    gy = g_ref * sin(tilt_ref)
  else
    gx = zero(L)
    gy = zero(L)
  end

  species = bunch.species
  tilde_m, _, beta_0 = BeamTracking.drift_params(species, p_over_q_ref)
  charge = chargeof(species)
  p0c = BeamTracking.R_to_pc(species, p_over_q_ref)
  mc2 = massof(species)
  n_steps, ds_step = BeamTracking.find_steps(tm, L)
  multipole_functions, multipole_parameters, multipole_normalized =
    runge_kutta_field(bmultipoleparams, L, p_over_q_ref)
  custom_functions, custom_parameters, custom_normalized =
    runge_kutta_custom_field(em_field_params)
  field_functions = (multipole_functions..., custom_functions...)
  field_parameters = (multipole_parameters..., custom_parameters...)
  field_normalized = (multipole_normalized..., custom_normalized...)
  length(field_functions) == length(field_parameters) == length(field_normalized) ||
    throw(DimensionMismatch("field functions, parameters, and normalization flags must have equal lengths"))

  # Time-dependent values in params are evaluated once, at the particle's
  # element-entrance time, by the common kernel path. They stay fixed during
  # all RK substeps.
  params = (beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps,
            gx, gy, field_functions, field_parameters, field_normalized, edge_params, edge_ksol, fin, fout)
  return push(kc, make_kernel_call(BeamTracking.rk4_kernel!, params))
end
