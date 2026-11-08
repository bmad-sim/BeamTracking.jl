@generated function c_light(::Type{T}) where {T}
  c = C_LIGHT
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    c = TS(c)
  end
  return :($c)
end

@generated function one_half(::Type{T}) where {T}
  h = 1/2
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    h = TS(h)
  end
  return :($h)
end

function float_type(::Type{T}) where {T}
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    return TS
  else # Includes TPSA, Duals, etc.
    return Float64
  end
end