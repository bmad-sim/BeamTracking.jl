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
const emptySA = SVector{1, Int32}[]

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

function rf_omega_calc(rfparams, circumference, species, p_over_q_ref)
  if rfparams.harmon_master
    tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(species, p_over_q_ref)
    return 2*pi*rfparams.harmon*C_LIGHT*beta_0/circumference
  else
    return 2*pi*rfparams.rf_frequency
  end
end

#---------------------------------------------------------------------------------------------------

function rf_phi0_calc(rfparams, species)
  chargeof(species) > 0 ? dphi = 0 : dphi = pi

  if rfparams.zero_phase == PhaseReference.BelowTransition
    return rfparams.phi0 + 0.5*pi + dphi
  elseif rfparams.zero_phase == PhaseReference.AboveTransition
    return rfparams.phi0 - 0.5*pi + dphi
  elseif rfparams.zero_phase == PhaseReference.Accelerating
    return rfparams.phi0 + dphi
  else
    error("RF parameter zero_phase value not set correctly.")
  end
end

#---------------------------------------------------------------------------------------------------

function rf_phi0_calc_old(rfparams, species)
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
    rf_step_calc(n_cell, L_active, rf_omega, L) -> n_cell_out, L_active_out

For an element of length `L`, calculate the number of RF cells (kicks) `n_cell_out` and the active length 
`L_active_out` given the input number of cells `n_cell` and the active length length `L_active`.

If `L_active` is negative, `L_active_out` is set to the element length `L`.
If `n_cell` is negative, `n_cell_out` is set so that the cell length is near half a wavelength.
If `L_active` is zero, `n_cell_out` is set to zero.
"""
function rf_step_calc(n_cell, L_active, rf_omega, L)  
  L_active < 0 ? L_act = L : L_act = L_active

  if n_cell < 0
    return round(rf_omega * L_act / (pi * C_LIGHT)), L_act
  elseif L_active == 0
    return 0, L_active
  else
    return n_cell, L_act
  end
end

#---------------------------------------------------------------------------------------------------

function bunch_dt_ref(tm, bunch, rfP, beamlineP, L)
  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)
   L / beta_gamma_to_v(beta_gamma_ref)
end


function bunch_dt_ref(tm::SaganCavity, bunch, rfP, beamlineP, L)
  if L == 0; return 0; end

  species = bunch.species
  p1_over_q_ref = beamlineP.beamline.p_over_q_ref
  rf_omega = rf_omega_calc(rfP, beamlineP.beamline.line[end].s_downstream, species, p1_over_q_ref)
  n_cell, L_active = rf_step_calc(tm.n_cell, tm.L_active, rf_omega, L)
  L_outer = (L - L_active) / 2
  E1_ref = R_to_E(species, p1_over_q_ref)
  dE_ref = beamlineP.dE_ref
  E0_ref = E1_ref - dE_ref
  dt_ref = L_outer/E_to_c_beta(species, E0_ref) + L_outer/E_to_c_beta(species, E1_ref)

  if n_cell == 0
    L_inner = L_active / 2
    dt_ref += L_inner/E_to_c_beta(species, E0_ref) + L_inner/E_to_c_beta(species, E1_ref)
  else
    for i_step = 1:n_cell
      E_now_ref = E0_ref + (i_step-0.5) * dE_ref / n_cell
      dt_ref += L_active / (n_cell * E_to_c_beta(species, E_now_ref))
    end
  end

  return dt_ref
end
