@generated function c_light(::Type{T}) where {T}
  c = C_LIGHT
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    c = TS(c)
  end
  return :($c)
end

function float_type(::Type{T}) where {T}
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    return TS
  else # Includes TPSA, Duals, etc.
    return Float64
  end
end