@generated function c_light(::Type{T}) where {T}
  c = C_LIGHT
  if T == Float32 || T == Float16
    c = T(c)
  end
  return :($c)
end