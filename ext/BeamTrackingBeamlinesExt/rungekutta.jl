# RungeKutta only constructs the body field integration kernel.
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

  if isactive(bendparams)
    (bendparams.edge1_int == 0 && bendparams.edge2_int == 0) ||
      error("edge1_int and edge2_int not yet handled for tracking")
    (bendparams.e1 == 0 && bendparams.e2 == 0) ||
      error("RungeKutta tracking does not support nonzero bend edge angles e1 or e2 because fringe tracking is not implemented")
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
            gx, gy, field_functions, field_parameters, field_normalized)
  return push(kc, make_kernel_call(BeamTracking.rk4_kernel!, params))
end
