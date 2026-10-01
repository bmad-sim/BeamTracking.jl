# Field values are flat tuples (Ex, Ey, Ez, Bx, By, Bz).
@inline add_fields(a, b) = map(+, a, b)
@inline scale_field(field, factor) = map(v -> v * factor, field)

"""
    zero_field(x, y, s, t, parameters=nothing)

Return `(Ex, Ey, Ez, Bx, By, Bz)` with zero values of the coordinate's numeric type.
"""
@inline function zero_field(x, y, s, t, parameters=nothing)
  v = zero(x)
  return (v, v, v, v, v, v)
end

"""
    multipole_field(x, y, s, t, parameters)

Evaluate magnetic multipoles from `(orders, normal, skew)` coefficients.
Element tracking obtains these coefficients exclusively from `BMultipoleParams`
during unpacking. Return `(Ex, Ey, Ez, Bx, By, Bz)` in the coefficients' units.
"""
@inline function multipole_field(x, y, s, t, parameters)
  orders, normal, skew = parameters
  bx, by = normalized_field(orders, normal, skew, x, y, 0)
  zero_field = zero(bx)
  bz = vifelse(orders[1] == 0, normal[1], zero_field)
  return (zero_field, zero_field, zero_field, bx, by, bz)
end

"""
    normalized_field_at(field_function, parameters, normalized, x, y, s, t, inv_rigidity)

Evaluate `field_function(x, y, s, t, parameters)`. `normalized` is `Val(true)` for an evaluator that
already returns normalized fields, and `Val(false)` for physical fields.

Parallel tuples of functions, parameters, and flags add contributions after
normalizing each independently. Empty tuples produce zero field.
"""
@inline function normalized_field_at(field_function, parameters, ::Val{normalized},
                                     x, y, s, t, inv_rigidity) where {normalized}
  field = field_function(x, y, s, t, parameters)
  if normalized
    return field
  else
    return scale_field(field, inv_rigidity)
  end
end

@inline function normalized_field_at(field_functions::Tuple, parameters::Tuple, normalized::Tuple,
                                     x, y, s, t, inv_rigidity)
  return _normalized_field_sum(field_functions, parameters, normalized, x, y, s, t, inv_rigidity)
end

@inline _normalized_field_sum(::Tuple{}, ::Tuple{}, ::Tuple{}, x, y, s, t, inv_rigidity) =
  zero_field(x, y, s, t)

@inline function _normalized_field_sum(field_functions::Tuple, parameters::Tuple, normalized::Tuple,
                                       x, y, s, t, inv_rigidity)
  field = normalized_field_at(first(field_functions), first(parameters), first(normalized),
                              x, y, s, t, inv_rigidity)
  # Avoid adding a zero field to the last contribution, especially for TPSA.
  length(field_functions) == 1 && return field
  return add_fields(field, _normalized_field_sum(Base.tail(field_functions), Base.tail(parameters),
                                                Base.tail(normalized), x, y, s, t, inv_rigidity))
end
