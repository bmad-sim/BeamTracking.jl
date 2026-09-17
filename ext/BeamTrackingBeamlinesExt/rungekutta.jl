# RungeKutta only constructs the body field integration kernel.
# Unpacking, reference-ramp, alignment, aperture, and callback
# are handled by the shared unpacking step.

@inline _unpack_field_parameter(value::NamedTuple{names}, context) where {names} =
  NamedTuple{names}(_unpack_field_parameter(Tuple(value), context))
@inline _unpack_field_parameter(value::Tuple, context) =
  map(item -> _unpack_field_parameter(item, context), value)
@inline _unpack_field_parameter(value::StaticArray, context) =
  map(item -> _unpack_field_parameter(item, context), value)
@inline _unpack_field_parameter(value, context) = deval(value, context)

@inline _scalarize_field_parameter(value::NamedTuple{names}) where {names} =
  NamedTuple{names}(_scalarize_field_parameter(Tuple(value)))
@inline _scalarize_field_parameter(value::Tuple) =
  map(_scalarize_field_parameter, value)
@inline _scalarize_field_parameter(value::StaticArray) =
  map(_scalarize_field_parameter, value)
@inline _scalarize_field_parameter(value) = scalarize(value)

@inline Beamlines.deval(field_source::MultipoleField, context) =
  BeamTracking._rebuild_multipole_field(
    field_source,
    _unpack_field_parameter(field_source.normal, context),
    _unpack_field_parameter(field_source.skew, context),
  )
@inline Beamlines.deval(field_source::FunctionalField, context) =
  BeamTracking._rebuild_functional_field(
    field_source,
    _unpack_field_parameter(field_source.parameters, context),
  )
@inline Beamlines.deval(field_source::SumField, context) =
  BeamTracking._rebuild_sum_field(
    field_source,
    map(item -> deval(item, context), field_source.field_sources),
  )

@inline Beamlines.scalarize(field_source::MultipoleField) =
  BeamTracking._rebuild_multipole_field(
    field_source,
    _scalarize_field_parameter(field_source.normal),
    _scalarize_field_parameter(field_source.skew),
  )
@inline Beamlines.scalarize(field_source::FunctionalField) =
  BeamTracking._rebuild_functional_field(
    field_source,
    _scalarize_field_parameter(field_source.parameters),
  )
@inline Beamlines.scalarize(field_source::SumField) =
  BeamTracking._rebuild_sum_field(
    field_source,
    map(scalarize, field_source.field_sources),
  )

@inline function runge_kutta_field(bmultipoleparams, L, p_over_q_ref)
  if !isactive(bmultipoleparams)
    return ZeroField()
  end

  mm = getfield(bmultipoleparams, :order)
  bn, bs = get_strengths(bmultipoleparams, L, p_over_q_ref)
  if mm isa Integer
    return MultipoleField(SA[mm], SA[bn], SA[bs]; normalized=true)
  end
  return MultipoleField(mm, bn, bs; normalized=true)
end

@inline configured_runge_kutta_field(::Nothing, element_field_source) = element_field_source

@inline function configured_runge_kutta_field(field_source_params::FieldSourceParams, element_field_source)
  field_source = field_source_params.field_source
  additional_field = field_source_params.additional_field
  if !isnothing(field_source) && !isnothing(additional_field)
    error("FieldSourceParams accepts either field_source or additional_field")
  elseif !isnothing(field_source)
    return field_source
  elseif !isnothing(additional_field)
    return SumField(element_field_source, additional_field)
  end
  return element_field_source
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
  field_source_params,
  L,
)
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
  element_field_source = runge_kutta_field(bmultipoleparams, L, p_over_q_ref)
  field_source = configured_runge_kutta_field(field_source_params, element_field_source)

  # Time-dependent values in params are evaluated once, at the particle's
  # element-entrance time, by the common kernel path. They stay fixed during
  # all RK substeps.
  params = (beta_0, tilde_m, charge, p0c, mc2, L, ds_step, n_steps,
            gx, gy, field_source)
  return push(kc, make_kernel_call(BeamTracking.rk4_kernel!, params))
end
