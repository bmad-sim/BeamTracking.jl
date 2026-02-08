Base.promote_rule(::Type{DefExpr{T}}, ::Type{TimeDependentParam}) where {T} = DefExpr{TimeDependentParam}

Beamlines.DefExpr{T}(a::TimeDependentParam) where {T} = DefExpr{T}(()->convert(T,a))
#DefExpr{T}(a::DefExpr) where {T} = DefExpr{T}(()->convert(T,a()))

function check_bl_bunch!(bl::Beamline, bunch::Bunch, notify::Bool=true)
  ref = getfield(bl, :ref)
  species_ref = getfield(bl, :species_ref)
  check_species!(species_ref, bunch, notify)
  check_p_over_q_ref!(bl, ref, bunch, notify)
  return
end

function check_species!(species_ref::Species, bunch::Bunch, notify=true)
  if isnullspecies(bunch.species)
    if isnullspecies(species_ref)
      error("Bunch species has not been set")
    else
      if notify
        println("Setting bunch.species = $species_ref (reference species from the Beamline)")
      end
      setfield!(bunch, :species, species_ref)
    end
  elseif !isnullspecies(species_ref) && species_ref != bunch.species && notify
    println("WARNING: The species of the bunch does NOT equal the reference species of the Beamline.")
  end
  return
end

function check_p_over_q_ref!(bl::Beamline, ref, bunch::Bunch, notify=true)
  t_ref = bunch.t_ref
  if isnan(bunch.p_over_q_ref)
    if isnothing(ref)
      if notify
        println("WARNING: Both the bunch and beamline do not have any set reference energy. If any LineElements have unnormalized fields stored as independent variables, there will be an error.")
      end
    else
      if bl isa Beamline
        p_over_q_ref = bl.p_over_q_ref
      else
        p_over_q_ref = ref
      end
      if notify
        if ref isa TimeDependentParam
          println("Setting bunch.p_over_q_ref = $(p_over_q_ref(t_ref)) (reference p_over_q_ref from the Beamline at t_ref = $t_ref)")
        else
          println("Setting bunch.p_over_q_ref = $p_over_q_ref (reference p_over_q_ref from the Beamline)")
        end
      end
      if p_over_q_ref isa TimeDependentParam
        setfield!(bunch, :p_over_q_ref, typeof(bunch.p_over_q_ref)(p_over_q_ref(t_ref)))
      else
        setfield!(bunch, :p_over_q_ref, typeof(bunch.p_over_q_ref)(p_over_q_ref))
      end
    end
  elseif !isnothing(ref)  && !(bl.p_over_q_ref ≈ bunch.p_over_q_ref) && !(bl.p_over_q_ref isa TimeDependentParam) && notify
    println("WARNING:The reference energy of the bunch does NOT equal the reference energy of the Beamline. 
              Normalized field strengths in tracking ALWAYS use the reference energy of the bunch.")
  end
  return
end

get_n_multipoles(::BMultipoleParams{T,N}) where {T,N} = N

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
@inline function get_strengths(bm, L, p_over_q_ref)
  if isconcretetype(eltype(bm.n))
    T = promote_type(eltype(bm.n),
                    typeof(L), typeof(p_over_q_ref)
    )
  else
    if bm.n isa AbstractArray
      T = promote_type(reduce(promote_type, typeof.(bm.n)), 
                      reduce(promote_type, typeof.(bm.s)),
                      reduce(promote_type, typeof.(bm.tilt)),
                      typeof(L), typeof(p_over_q_ref)
      )
    else
      T = promote_type(typeof(bm.n), 
                      typeof(bm.s),
                      typeof(bm.tilt),
                      typeof(L), typeof(p_over_q_ref)
      )
    end
  end
  n = T.(make_static(bm.n))
  s = T.(make_static(bm.s))
  tilt = T.(make_static(bm.tilt))
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  np = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/p_over_q_ref, np) 
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp) 
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

@inline function get_integrated_strengths(bm, L, p_over_q_ref)
  if isconcretetype(eltype(bm.n))
    T = promote_type(eltype(bm.n),
                    typeof(L), typeof(p_over_q_ref)
    )
  else
    if bm.n isa AbstractArray
      T = promote_type(reduce(promote_type, typeof.(bm.n)), 
                      reduce(promote_type, typeof.(bm.s)),
                      reduce(promote_type, typeof.(bm.tilt)),
                      typeof(L), typeof(p_over_q_ref)
      )
    else
      T = promote_type(typeof(bm.n), 
                      typeof(bm.s),
                      typeof(bm.tilt),
                      typeof(L), typeof(p_over_q_ref)
      )
    end
  end
  n = T.(make_static(bm.n))
  s = T.(make_static(bm.s))
  tilt = T.(make_static(bm.tilt))
  order = bm.order
  normalized = bm.normalized
  integrated = bm.integrated
  np = @. n*cos(order*tilt) + s*sin(order*tilt)
  sp = @. -n*sin(order*tilt) + s*cos(order*tilt)
  np = @. ifelse(!normalized, np/p_over_q_ref, np) 
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp) 
  np = @. ifelse(!integrated, np*L, np)
  sp = @. ifelse(!integrated, sp*L, sp)
  return np, sp
end

#---------------------------------------------------------------------------------------------------

function rf_omega(rfparams, circumference, species, p_over_q_ref)
  if rfparams.harmon_master
    tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(species, p_over_q_ref)
    return 2*pi*rfparams.harmon*C_LIGHT*beta_0/circumference
  else
    return 2*pi*rfparams.rf_frequency
  end
end

#---------------------------------------------------------------------------------------------------

function rf_phi0(rfparams)
  if rfparams.zero_phase == PhaseReference.BelowTransition
    return rfparams.phi0 + 0.5*pi
  elseif rfparams.zero_phase == PhaseReference.AboveTransition
    return rfparams.phi0 - 0.5*pi
  elseif rfparams.zero_phase == PhaseReference.Accelerating
    return rfparams.phi0
  else
    error("RF parameter zero_phase value not set correctly.")
  end
end