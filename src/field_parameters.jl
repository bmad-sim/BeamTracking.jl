"""
    field_parameter_leaf(::Type) -> Bool

Choose whether field-source preparation treats a type as an opaque leaf.
Ordinary structs are traversed through their fields in declaration order.
Arrays other than `SVector`, numbers, functions, and dynamic parameter wrappers
are leaves; operations handle them without inspecting their internal storage.
"""
@inline field_parameter_leaf(::Type{T}) where {T} = isprimitivetype(T) || fieldcount(T) == 0

@inline field_parameter_leaf(::Type{<:Union{Number,AbstractArray,Function,Type,AbstractString,
                                          Symbol,TimeDependentParam,TimeFunction,BatchParam,
                                          _LoweredBatchParam,SIMD.Vec}}) = true
@inline field_parameter_leaf(::Type{<:SVector}) = false

"""
    rebuild_field_source(original, children::Tuple)

Reconstruct a field source or parameter container after transforming its fields.
The default calls the unparameterized constructor with the children in field
order, allowing parameter types to change. Override
this method for constructors that do not follow that convention. Mutable objects
are reconstructed, not updated in place; cyclic object graphs are not supported.
"""
@generated function rebuild_field_source(original::T, children::Tuple) where {T}
  constructor = Base.typename(T).wrapper
  return :($constructor(children...))
end

@inline rebuild_field_source(::NamedTuple{names}, children::Tuple) where {names} =
  NamedTuple{names}(children)
@inline rebuild_field_source(::Tuple, children::Tuple) = children
@inline rebuild_field_source(::SVector{N}, children::Tuple) where {N} = SVector{N}(children)
@inline function rebuild_field_source(source::MultipoleField, children::Tuple)
  orders, normal, skew = children
  return MultipoleField{typeof(source.orders),typeof(normal),typeof(skew)}(
    orders, normal, skew,
  )
end
# Preserve the already flattened sum without repeating constructor normalization.
@inline rebuild_field_source(::SumField, children::Tuple) =
  SumField{typeof(only(children))}(only(children))

# Emit direct field accesses so inference follows nested heterogeneous structs
# without an intermediate recursive traversal of reflected field tuples.
@generated function map_field_parameters(operation::F, value::T, args::Vararg{Any,N}) where {F,T,N}
  isprimitivetype(T) && return :(operation(value, args...))
  children = [:(map_field_parameters(operation, getfield(value, $i), args...))
              for i in 1:fieldcount(T)]
  return quote
    $(Expr(:meta, :inline))
    if field_parameter_leaf(T)
      operation(value, args...)
    else
      rebuild_field_source(value, tuple($(children...)))
    end
  end
end

@generated function any_field_parameter(predicate, value::T) where {T}
  isprimitivetype(T) && return :(predicate(value))
  result = :(false)
  for i in fieldcount(T):-1:1
    result = :(any_field_parameter(predicate, getfield(value, $i)) || $result)
  end
  return quote
    $(Expr(:meta, :inline))
    field_parameter_leaf(T) ? predicate(value) : $result
  end
end

@inline map_field_parameters(operation, source::FunctionalField, args...) =
  FunctionalField(source.evaluator, map_field_parameters(operation, source.parameters, args...))
@inline any_field_parameter(predicate, source::FunctionalField) =
  any_field_parameter(predicate, source.parameters)
@inline map_field_parameters(operation, values::SVector, args...) =
  rebuild_field_source(values, map_field_parameters(operation, Tuple(values), args...))
@inline any_field_parameter(predicate, values::SVector) =
  any_field_parameter(predicate, Tuple(values))

# Dynamic static vectors become tuples before particle evaluation, so different
# closure types and SIMD lanes do not have to share one array element type.
@inline function map_field_parameters(::typeof(batch_lower), values::SVector)
  children = map_field_parameters(batch_lower, Tuple(values))
  return any_field_parameter(_is_batch_parameter, values) ? children :
         rebuild_field_source(values, children)
end
@inline function map_field_parameters(::typeof(time_lower), values::SVector)
  children = map_field_parameters(time_lower, Tuple(values))
  return any_field_parameter(_is_time_parameter, values) ? children :
         rebuild_field_source(values, children)
end
@inline _is_batch_parameter(value) = value isa BatchParam
@inline _is_time_parameter(value) = value isa TimeDependentParam

# Only field-source kernel arguments enter this traversal. Other kernel
# parameters retain their existing batch/time semantics.
struct _PreparedField{S}
  source::S
end
@inline (field::_PreparedField)(x, y, z, s) = field.source(x, y, z, s)
const _TraversedField = Union{MultipoleField,FunctionalField,SumField,_PreparedField}

@inline batch_lower(source::_TraversedField) = map_field_parameters(batch_lower, source)
@inline time_lower(source::_TraversedField) = map_field_parameters(time_lower, source)
@inline beval(source::_TraversedField, i) = map_field_parameters(beval, source, i)
@inline teval(source::_TraversedField, t) = map_field_parameters(teval, source, t)
@inline static_batchcheck(source::_TraversedField) = any_field_parameter(static_batchcheck, source)
@inline static_timecheck(source::_TraversedField) = any_field_parameter(static_timecheck, source)

Adapt.@adapt_structure _PreparedField
