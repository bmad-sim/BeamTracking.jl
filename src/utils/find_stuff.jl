# Temporary disgusting solution for callbacks - Symplectic
function _tilde_m_position(f)
  for m in methods(f)
    i = findfirst(==(:tilde_m), Base.method_argnames(m))
    i === nothing || return i - 1 - 3      # skip #self#
  end
  return 0
end

find_tilde_m(::Type{K}) where {K} = 0     # fallback

for k in (implicit_integrator!, exact_drift!, sks_multipole!,
          dkd_multipole!, bkb_multipole!, mkm_quadrupole!)
  @eval find_tilde_m(::Type{typeof($k)}) = $(_tilde_m_position(k))
end

@generated function beta_gamma_to_v(beta_gamma::T) where {T}
  clight = C_LIGHT
  TS = ForwardDiff.valtype(T)
  if TS == Float32 || TS == Float16
    clight = TS(C_LIGHT)
  end
  return :($clight*beta_gamma/sqrt(1+beta_gamma^2))
end

@generated function compute_dt_ref(s, ker::K, params) where {K}
  idx = find_tilde_m(ker)
  if idx == 0
    return quote
      return 0
    end
  else
    return quote
      beta_gamma_ref = 1 / params[$idx]
      return s / beta_gamma_to_v(beta_gamma_ref)
    end
  end
end
