function check_bl_bunch!(bl::Beamline, bunch::Bunch, notify::Bool=true)
  R_ref = getfield(bl, :R_ref)
  species_ref = getfield(bl, :species_ref)
  check_species!(species_ref, bunch, notify)
  check_R_ref!(R_ref, bunch, notify)
  return
end

function check_species!(species_ref::Species, bunch::Bunch, notify=true)
  if isnullspecies(bunch.species)
    if isnullspecies(species_ref)
      error("Bunch species has not been set")
    else
      if notify
        @info "Setting bunch.species = $species_ref (reference species from the Beamline)"
      end
      setfield!(bunch, :species, species_ref)
    end
  elseif !isnullspecies(species_ref) && species_ref != bunch.species && notify
    @warn "The species of the bunch does NOT equal the reference species of the Beamline."
  end
  return
end

function check_R_ref!(R_ref, bunch::Bunch, notify=true)
  if isnan(bunch.R_ref)
    if isnothing(R_ref)
      if notify
        @warn "Both the bunch and beamline do not have any set R_ref. If any LineElements have unnormalized fields stored as independent variables, there will be an error."
      end
    else
      if notify
        @info "Setting bunch.R_ref = $R_ref (reference R_ref from the Beamline)"
      end
      setfield!(bunch, :R_ref, typeof(bunch.R_ref)(R_ref)) #R_ref = R_ref
    end
  elseif !isnothing(R_ref)  && !(R_ref ≈ bunch.R_ref) && notify
    @warn "The reference energy of the bunch does NOT equal the reference energy of the Beamline. 
    Normalized field strengths in tracking ALWAYS use the reference energy of the bunch."
  end
  return
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
@inline function get_strengths(bm, L, R_ref)
  n = make_static(bm.n)
  s = make_static(bm.s)
  tilt = bm.tilt
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  # Make all the same type
  T = promote_type(Base.promote_op(/, typeof(n), typeof(R_ref)), 
                   Base.promote_op(/, typeof(n), typeof(L)))
  np::T = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp::T = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/R_ref, np) 
  sp = @. ifelse(!normalized, sp/R_ref, sp) 
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

@inline function get_integrated_strengths(bm, L, R_ref)
  n = make_static(bm.n)
  s = make_static(bm.s)
  tilt = bm.tilt
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  # Make all the same type
  T = promote_type(Base.promote_op(/, typeof(n), typeof(R_ref)), 
                   Base.promote_op(/, typeof(n), typeof(L)))
  np::T = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp::T = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/R_ref, np) 
  sp = @. ifelse(!normalized, sp/R_ref, sp) 
  np = @. ifelse(!integrated, np*L, np)
  sp = @. ifelse(!integrated, sp*L, sp)
  return np, sp
end


