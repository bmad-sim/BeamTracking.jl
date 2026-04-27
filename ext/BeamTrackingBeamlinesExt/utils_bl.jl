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

#---------------------------------------------------------------------------------------------------

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

#---------------------------------------------------------------------------------------------------

get_n_multipoles(::BMultipoleParams{T,N}) where {T,N} = N

make_static(a::StaticArray) = SVector(a)
make_static(a) = a

#---------------------------------------------------------------------------------------------------

"""
    get_strengths(bm, L, p_over_q_ref) -> Kn, Ks
Get non-integrated magnetic multipole strength arrays. Also see get_integrated_strengths.

(Kn' + im*Ks') = (Kn + im*Ks)*exp(-im*order*tilt)

# Rotation:
Kn' = Kn*cos(order*tilt) + Ks*sin(order*tilt)
Ks' = Kn*-sin(order*tilt) + Ks*cos(order*tilt)

This works for both BMultipole and BMultipoleParams. Branchless bc SIMD -> basically 
no loss in computing both but benefit of branchless.
"""
@inline function get_strengths(bm, L, p_over_q_ref)
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
  np = @. ifelse(!normalized, np/p_over_q_ref, np) 
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp) 
  np = @. ifelse(integrated, np/L, np)
  sp = @. ifelse(integrated, sp/L, sp)
  return np, sp
end

@inline function get_integrated_strengths(bm, L, p_over_q_ref)
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
  np = @. ifelse(!normalized, np/p_over_q_ref, np) 
  sp = @. ifelse(!normalized, sp/p_over_q_ref, sp) 
  np = @. ifelse(!integrated, np*L, np)
  sp = @. ifelse(!integrated, sp*L, sp)
  return np, sp
end

#---------------------------------------------------------------------------------------------------

function rf_omega_calc(rfparams, beamlineparams)
  if rfparams.rate_meaning == RateMeaning.Harmon
    beamline = beamlineparams.beamline
    p_over_q_ref = beamline.p_over_q_ref
    species = beamline.species_ref
    circumference = beamline.line[end].s_downstream
    v = R_to_v(species, p_over_q_ref)
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

function fringe_in(f::Fringe.T)
  if f == Fringe.BothEnds || f == Fringe.EntranceEnd
    return Val{true}()
  else
    return Val{false}()
  end
end

function fringe_out(f::Fringe.T)
  if f == Fringe.BothEnds || f == Fringe.ExitEnd
    return Val{true}()
  else
    return Val{false}()
  end
end

#---------------------------------------------------------------------------------------------------
"""
    rf_step_calc(num_cells, L_active, rf_omega, L) -> num_cells_out, L_active_out

For an element of length `L`, calculate the number of RF cells (kicks) `num_cells_out` and the active length 
`L_active_out` given the input number of cells `num_cells` and the active length length `L_active`.

If `L_active` is negative, `L_active_out` is set to the element length `L`.
If `num_cells` is negative, `num_cells_out` is set so that the cell length is near half a wavelength.
If `L_active` is zero, `num_cells_out` is set to zero.
"""
function rf_step_calc(num_cells, L_active, rf_omega, L)  
  L_active < 0 ? L_act = L : L_act = L_active

  if num_cells < 0
    return round(rf_omega * L_act / (pi * C_LIGHT)), L_act
  elseif L_active == 0
    return 0, L_active
  else
    return num_cells, L_act
  end
end
