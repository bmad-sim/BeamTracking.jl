function check_Brho(Brho_ref, bunch::Bunch)
  if isnan(bunch.Brho_ref) && isnan(Brho_ref)
    @warn "Both the bunch and beamline do not have any set Brho_ref. If any LineElements have unnormalized fields stored as independent variables, tracking results will be NaNs"
    setfield!(bunch, :Brho_ref, typeof(bunch.Brho_ref)(Brho_ref))
    Bρ = (Brho_ref, Brho_ref)
  else
    Bρ = (bunch.Brho_ref, Brho_ref)
    setfield!(bunch, :Brho_ref, typeof(bunch.Brho_ref)(Brho_ref))
  end
  return Bρ
end

get_n_multipoles(::BMultipoleParams{T,N}) where {T,N} = N
function get_n_multipoles(b::Beamlines.BitsBMultipoleParams{T,N}) where {T,N}
  n = 0
  i = 1
  while i <= N && b.order[i] >= 0
    n += 1
    i += 1
  end
  return n
end

make_static(a::StaticArray) = SVector(a)
make_static(a) = a

"""

(Kn' + im*Ks') = (Kn + im*Ks)*exp(-im*order*tilt)

# Rotation:
Kn' = Kn*cos(order*tilt) + Ks*sin(order*tilt)
Ks' = Kn*-sin(order*tilt) + Ks*cos(order*tilt)

This works for both BMultipole and BMultipoleParams. Branchless bc SIMD -> basically 
no loss in computing both but benefit of branchless.
"""
@inline function get_strengths(bm, L, Brho_ref)
  n = make_static(bm.n)
  s = make_static(bm.s)
  tilt = bm.tilt
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  # Make all the same type
  T = promote_type(Base.promote_op(/, typeof(n), typeof(Brho_ref)), 
                   Base.promote_op(/, typeof(n), typeof(L)))
  np::T = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp::T = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/Brho_ref, np) 
  sp = @. ifelse(!normalized, sp/Brho_ref, sp) 
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

@inline function get_integrated_strengths(bm, L, Brho_ref)
  n = make_static(bm.n)
  s = make_static(bm.s)
  tilt = bm.tilt
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  # Make all the same type
  T = promote_type(Base.promote_op(/, typeof(n), typeof(Brho_ref)), 
                   Base.promote_op(/, typeof(n), typeof(L)))
  np::T = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp::T = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/Brho_ref, np) 
  sp = @. ifelse(!normalized, sp/Brho_ref, sp) 
  np = @. ifelse(!integrated, np*L, np)
  sp = @. ifelse(!integrated, sp*L, sp)
  return np, sp
end


