@inline function _rebuild_multipole_field(field_source::MultipoleField{M,KN,KS,N}, normal, skew) where {M,KN,KS,N}
  return MultipoleField{typeof(field_source.orders),typeof(normal),typeof(skew),N}(
    field_source.orders,
    normal,
    skew,
  )
end

@inline batch_lower(field_source::MultipoleField) =
  _rebuild_multipole_field(field_source, batch_lower(field_source.normal), batch_lower(field_source.skew))
@inline time_lower(field_source::MultipoleField) =
  _rebuild_multipole_field(field_source, time_lower(field_source.normal), time_lower(field_source.skew))
@inline num_lower(::Type{T}, field_source::MultipoleField) where {T<:Union{Float32,Float16}} =
  _rebuild_multipole_field(field_source, num_lower(T, field_source.normal), num_lower(T, field_source.skew))
@inline static_batchcheck(field_source::MultipoleField) =
  static_batchcheck(field_source.normal) || static_batchcheck(field_source.skew)
@inline static_timecheck(field_source::MultipoleField) =
  static_timecheck(field_source.normal) || static_timecheck(field_source.skew)
@inline beval(field_source::MultipoleField, i) =
  _rebuild_multipole_field(field_source, beval(field_source.normal, i), beval(field_source.skew, i))
@inline teval(field_source::MultipoleField, t) =
  _rebuild_multipole_field(field_source, teval(field_source.normal, t), teval(field_source.skew, t))

@inline _rebuild_functional_field(field_source::FunctionalField{F,P,N}, parameters) where {F,P,N} =
  FunctionalField(field_source.evaluator, parameters; normalized=N)

@inline batch_lower(field_source::FunctionalField) =
  _rebuild_functional_field(field_source, batch_lower(field_source.parameters))
@inline time_lower(field_source::FunctionalField) =
  _rebuild_functional_field(field_source, time_lower(field_source.parameters))
@inline num_lower(::Type{T}, field_source::FunctionalField) where {T<:Union{Float32,Float16}} =
  _rebuild_functional_field(field_source, num_lower(T, field_source.parameters))
@inline static_batchcheck(field_source::FunctionalField) =
  static_batchcheck(field_source.parameters)
@inline static_timecheck(field_source::FunctionalField) =
  static_timecheck(field_source.parameters)
@inline beval(field_source::FunctionalField, i) =
  _rebuild_functional_field(field_source, beval(field_source.parameters, i))
@inline teval(field_source::FunctionalField, t) =
  _rebuild_functional_field(field_source, teval(field_source.parameters, t))

@inline _rebuild_sum_field(field_source::SumField, field_sources) =
  SumField{typeof(field_sources)}(field_sources)

@inline batch_lower(field_source::SumField) =
  _rebuild_sum_field(field_source, batch_lower(field_source.field_sources))
@inline time_lower(field_source::SumField) =
  _rebuild_sum_field(field_source, time_lower(field_source.field_sources))
@inline num_lower(::Type{T}, field_source::SumField) where {T<:Union{Float32,Float16}} =
  _rebuild_sum_field(field_source, num_lower(T, field_source.field_sources))
@inline static_batchcheck(field_source::SumField) =
  static_batchcheck(field_source.field_sources)
@inline static_timecheck(field_source::SumField) =
  static_timecheck(field_source.field_sources)
@inline beval(field_source::SumField, i) =
  _rebuild_sum_field(field_source, beval(field_source.field_sources, i))
@inline teval(field_source::SumField, t) =
  _rebuild_sum_field(field_source, teval(field_source.field_sources, t))
