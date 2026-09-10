@generated function c_light(::Type{T}) where {T}
  c = C_LIGHT
  if T == Float32 || T == Float16
    c = T(c)
  end
  return :($c)
end

function float_type(::Type{T}) where {T}
  if T == Float32 || T == Float16
    return T
  else # Includes TPSA, Duals, etc.
    return Float64
  end
end