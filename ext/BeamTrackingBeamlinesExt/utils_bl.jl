Base.promote_rule(::Type{DefExpr{T}}, ::Type{TimeDependentParam}) where {T} = DefExpr{TimeDependentParam}
Base.promote_rule(::Type{DefExpr{T}}, ::Type{BatchParam}) where {T} = DefExpr{promote_type(BatchParam, T)}
Beamlines.DefExpr{T}(a::TimeDependentParam) where {T} = DefExpr{T}(()->convert(T,a))
Beamlines.DefExpr{BatchParam}(a::BatchParam) = DefExpr{BatchParam}(()->a)

"""
    TrackReference

The references a `track!` call tracks against, resolved once by `resolve_reference!`:

- `species_ref`: species the Beamline fields and RF settings are designed for
- `p_over_q_ref`: reference p_over_q_ref the Bunch is tracked with in every element (Number or 
  TimeDependentParam). This is the Beamline p_over_q_ref times `field_scale` (same magnetic rigidity 
  magnitude, sign of the Bunch charge), or the Bunch's own if the Beamline has none set.
- `field_scale`: Beamline p_over_q_ref over the tracking `p_over_q_ref`: +1, or -1 if the Bunch 
  species has the opposite charge sign of `species_ref`
"""
struct TrackReference{P}
  species_ref::Species
  p_over_q_ref::P
  field_scale::Int
end

"""
    resolve_reference!(bunch::Bunch, species_ref, p_over_q_ref, notify::Bool=true) -> TrackReference

Resolves the reference species and energy between the Bunch and the Beamline. `species_ref` and 
`p_over_q_ref` are the Beamline's, or `nothing`/`GetError` if not set. If one of the Bunch or 
Beamline does not have a species or energy set, it is inherited from the other. If both have 
a reference energy set but they differ, the Bunch is ramped to the Beamline reference energy when 
entering the first element. If both have a species set but they differ, the fields of the Beamline 
are held fixed and the Bunch is tracked with the same magnitude of p_over_q_ref (magnetic rigidity) 
as the Beamline.
"""
function resolve_reference!(bunch::Bunch, species_ref, p_over_q_ref, notify::Bool=true)
  unset(x) = isnothing(x) || x isa Beamlines.GetError

  # Species
  if isnullspecies(bunch.species)
    unset(species_ref) && error("Neither the Bunch nor the Beamline has a species set")
    notify && println("Setting bunch.species = $species_ref (reference species from the Beamline)")
    setfield!(bunch, :species, species_ref)
  elseif unset(species_ref)
    species_ref = bunch.species
  elseif species_ref != bunch.species
    unset(p_over_q_ref) && error("The Bunch species ($(bunch.species)) differs from the Beamline reference species ($species_ref), but the Beamline has no reference energy set, so its normalized fields have no physical value")
    notify && println("Tracking $(bunch.species) through a Beamline with reference species $species_ref: the fields of the Beamline are held fixed")
  end
  field_scale = sign(chargeof(bunch.species)) == sign(chargeof(species_ref)) ? 1 : -1

  # Energy
  if unset(p_over_q_ref)
    isnan(bunch.p_over_q_ref) && error("Neither the Bunch nor the Beamline has a reference energy set")
    p_over_q_ref = bunch.p_over_q_ref
  else
    # Same magnetic rigidity magnitude as the Beamline, with the sign of the Bunch charge
    p_over_q_ref = field_scale == 1 ? p_over_q_ref : -p_over_q_ref
    if isnan(bunch.p_over_q_ref)
      t_ref = bunch.t_ref
      p0 = p_over_q_ref isa TimeDependentParam ? p_over_q_ref(t_ref) : p_over_q_ref
      notify && println("Setting bunch.p_over_q_ref = $p0 (reference p_over_q_ref from the Beamline at t_ref = $t_ref)")
      setfield!(bunch, :p_over_q_ref, typeof(bunch.p_over_q_ref)(p0))
    end
  end
  sign(bunch.p_over_q_ref) == sign(chargeof(bunch.species)) || error("bunch.p_over_q_ref = $(bunch.p_over_q_ref) has the wrong sign for species $(bunch.species)")
  # If both are set but differ, the Bunch is ramped to p_over_q_ref when entering the first element

  return TrackReference(species_ref, p_over_q_ref, field_scale)
end

#---------------------------------------------------------------------------------------------------

get_n_multipoles(::BMultipoleParams{T,N}) where {T,N} = N
get_n_multipoles(::EMultipoleParams{T,N}) where {T,N} = N

make_static(a::StaticArray) = SVector(a)
make_static(a::SizedArray) = SVector(a)
make_static(a) = a

#---------------------------------------------------------------------------------------------------

"""
    get_strengths(bm, L, p_over_q_ref, field_scale) -> Kn, Ks
Get non-integrated magnetic multipole strength arrays. Also see get_integrated_strengths.

`p_over_q_ref` is the tracking (Bunch) reference. Unnormalized strengths are divided by it, and
normalized strengths are multiplied by `field_scale` (see `TrackReference`), so that the physical 
fields are the same for any tracked species.

(Kn' + im*Ks') = (Kn + im*Ks)*exp(-im*order*tilt)

# Rotation:
Kn' = Kn*cos(order*tilt) + Ks*sin(order*tilt)
Ks' = Kn*-sin(order*tilt) + Ks*cos(order*tilt)

This works for both BMultipole and BMultipoleParams. Branchless bc SIMD -> basically 
no loss in computing both but benefit of branchless.
"""
@inline function get_strengths(bm, L, p_over_q_ref, field_scale)
  bmn = getfield(bm, :n)
  bms = getfield(bm, :s)
  bmtilt = getfield(bm, :tilt)
  if isconcretetype(eltype(bmn))
    T = promote_type(eltype(bmn),
                    typeof(L), typeof(p_over_q_ref)
    )
  else
    if bmn isa AbstractArray
      T = promote_type(reduce(promote_type, typeof.(bmn)), 
                      reduce(promote_type, typeof.(bms)),
                      reduce(promote_type, typeof.(bmtilt)),
                      typeof(L), typeof(p_over_q_ref)
      )
    else
      T = promote_type(typeof(bmn), 
                      typeof(bms),
                      typeof(bmtilt),
                      typeof(L), typeof(p_over_q_ref)
      )
    end
  end
  n = T.(make_static(bmn))
  s = T.(make_static(bms))
  tilt = T.(make_static(bmtilt))
  order = getfield(bm, :order)
  normalized = getfield(bm, :normalized)
  integrated = getfield(bm, :integrated)
  np = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/p_over_q_ref, np*field_scale)
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp*field_scale)
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

@inline function get_integrated_strengths(bm, L, p_over_q_ref, field_scale)
  bmn = getfield(bm, :n)
  bms = getfield(bm, :s)
  bmtilt = getfield(bm, :tilt)
  if isconcretetype(eltype(bmn))
    T = promote_type(eltype(bmn),
                    typeof(L), typeof(p_over_q_ref)
    )
  else
    if bmn isa AbstractArray
      T = promote_type(reduce(promote_type, typeof.(bmn)), 
                      reduce(promote_type, typeof.(bms)),
                      reduce(promote_type, typeof.(bmtilt)),
                      typeof(L), typeof(p_over_q_ref)
      )
    else
      T = promote_type(typeof(bmn), 
                      typeof(bms),
                      typeof(bmtilt),
                      typeof(L), typeof(p_over_q_ref)
      )
    end
  end
  n = T.(make_static(bmn))
  s = T.(make_static(bms))
  tilt = T.(make_static(bmtilt))
  order = getfield(bm, :order)
  normalized = getfield(bm, :normalized)
  integrated = getfield(bm, :integrated)
  np = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/p_over_q_ref, np*field_scale)
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp*field_scale)
  np = @. ifelse(!integrated, np*L, np)
  sp = @. ifelse(!integrated, sp*L, sp)
  return np, sp
end

@inline function get_e_strengths(em, L)
  emn = getfield(em, :n)
  ems = getfield(em, :s)
  emtilt = getfield(em, :etilt)
  if isconcretetype(eltype(emn))
    T = promote_type(eltype(emn),
                    typeof(L)
    )
  else
    if emn isa AbstractArray
      T = promote_type(reduce(promote_type, typeof.(emn)), 
                      reduce(promote_type, typeof.(ems)),
                      reduce(promote_type, typeof.(emtilt)),
                      typeof(L)
      )
    else
      T = promote_type(typeof(emn), 
                      typeof(ems),
                      typeof(emtilt),
                      typeof(L)
      )
    end
  end
  n = T.(make_static(emn))
  s = T.(make_static(ems))
  tilt = T.(make_static(emtilt))
  order = getfield(em, :order)
  integrated = getfield(em, :integrated)
  np = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

#---------------------------------------------------------------------------------------------------

function rf_omega_calc(rfparams, beamlineparams, ref)
  if rfparams.rate_meaning == RateMeaning.Harmon
    circumference = beamlineparams.beamline.line[end].s_downstream
    v = R_to_v(ref.species_ref, ref.p_over_q_ref) # Same |p_over_q_ref| as the Beamline, sign does not matter
    return 2*pi*rfparams.rate*v/circumference
  else
    return 2*pi*rfparams.rate
  end
end

#---------------------------------------------------------------------------------------------------

function rf_phi0_calc(rfparams, species)
  dphi = chargeof(species) > 0 ? 0 : pi

  if rfparams.zero_phase == PhaseRef.BelowTransition
    return rfparams.phi0 + pi/2 + dphi
  elseif rfparams.zero_phase == PhaseRef.AboveTransition
    return rfparams.phi0 - pi/2 + dphi
  elseif rfparams.zero_phase == PhaseRef.Accelerating
    return rfparams.phi0 + dphi
  else
    error("RF parameter zero_phase value not set correctly.")
  end
end

#---------------------------------------------------------------------------------------------------
"""
    rf_step_calc(n_cells, L_active, rf_omega, L) -> n_cells_out, L_active_out

For an element of length `L`, calculate the number of RF cells (kicks) `n_cells_out` and the active length 
`L_active_out` given the input number of cells `n_cells` and the active length length `L_active`.

If `L_active` is negative, `L_active_out` is set to the element length `L`.
If `n_cells` is negative, `n_cells_out` is set so that the cell length is near half a wavelength.
If `L_active` is zero, `n_cells_out` is set to zero.
"""
function rf_step_calc(n_cells, L_active, rf_omega, L)  
  L_active < 0 ? L_act = L : L_act = L_active

  if n_cells < 0
    return round(rf_omega * L_act / (pi * C_LIGHT)), L_act
  elseif L_active == 0
    return 0, L_active
  else
    return n_cells, L_act
  end
end

#---------------------------------------------------------------------------------------------------

function drift_params(species::Species, p_over_q_ref)
  beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
  tilde_m = 1/beta_gamma_0
  gamsqr_0 = @FastGTPSA 1+beta_gamma_0*beta_gamma_0
  beta_0 = @FastGTPSA beta_gamma_0/sqrt(gamsqr_0)
  return tilde_m, gamsqr_0, beta_0
end