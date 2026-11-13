#---------------------------------------------------------------------------------------------------
# universal! for SaganCavity tracking.
# its own universal! is needed because it does its own thing to the reference time
function universal!(
  bunch,
  tm::SaganCavity,
  fullele,
  universalparams,
  alignmentparams,
  bendparams,
  bmultipoleparams,
  patchparams,
  apertureparams,
  rfparams,
  beamlineparams,
  mapparams,
  fourpotentialparams,
  emultipoleparams,
  emfieldparams,
  ibsparams,
  ref,
  ::Val{ramp_particle_energy_without_rf}, 
  ::Val{ramp_update_each_particle},
  ::Val{shift_time},
  batch_start,
  sink;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle, shift_time}
  ele = LoweredElement(universalparams, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams, rfparams, 
                       beamlineparams, mapparams, fourpotentialparams, emultipoleparams, emfieldparams, ibsparams, ref)
  L = ele.universalparams.L
  coords = bunch.coords
  !isactive(mapparams) || error("SaganCavity Tracking through element $(fullele.name) with MapParams is undefined")
  !isactive(patchparams) || error("SaganCavity Tracking through element $(fullele.name) with PatchParams is undefined")
  !isactive(bendparams) || error("SaganCavity Tracking through element $(fullele.name) with BendParams is undefined")
  !isactive(fourpotentialparams) || error("SaganCavity Tracking through element $(fullele.name) with FourPotentialParams is undefined")
  !isactive(emultipoleparams) || error("SaganCavity Tracking through element $(fullele.name) with EMultipoleParams is undefined")
  isactive(rfparams) || error("SaganCavity Tracking through element $(fullele.name) without RFParams is undefined")
  universalparams.n_steps == 1 || error("SaganCavity is not compatible with elements with n_steps > 1 (its kicks are set by n_cells). Try n_steps = 1 instead")

  p_over_q_ref = ref.p_over_q_ref
  # A ramped reference is not supported by SaganCavity (see "Status" in AGENTS.md); when preparing it
  # would also need the reference momentum at the call
  !(sink isa Bunch) && p_over_q_ref isa TimeDependentParam &&
    error("SaganCavity with a time dependent reference energy is not supported by prepared tracking")
  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)

  # Compute information about reference coordinate system:
  t_enter = bunch.t_ref
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref(t_enter))
  else
    beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref)
  end
  g = isnothing(bendparams) ? (0,0) : reverse((bendparams.g_ref .* sincos(bendparams.tilt_ref)))
  ds_step = ele.universalparams.ds_step
  # Reference energy change. For a Bunch species other than the reference species, the Bunch 
  # reference energy changes between the energies of the Bunch species at the upstream and this 
  # Beamline's p_over_q_ref (the same |Bρ| as the design), the RF phase is relative to the design 
  # particle (reference species, see design_params), and z is re-referenced to it at the exit
  species = bunch.species
  p1_over_q_ref = p_over_q_ref
  E1_ref = R_to_E(species, p1_over_q_ref)
  if shift_time
    dE_ref = E1_ref - R_to_E(species, ele.ref.field_scale*p1_over_q_ref - beamlineparams.dp_over_q_ref)
    species_d = ele.ref.species_ref
    E1_d = R_to_E(species_d, ele.ref.field_scale*p1_over_q_ref)
    dE_d = beamlineparams.dE_ref
    design_params = (massof(species_d), sqrt((E1_d - dE_d)^2 - massof(species_d)^2), dE_d)
  else
    dE_ref = beamlineparams.dE_ref
    design_params = nothing
  end

  # reference time change
  if L != 0
    rf_omega = rf_omega_calc(rfparams, beamlineparams, ele.ref)
    n_cells, L_active = rf_step_calc(tm.n_cells, tm.L_active, rf_omega, L)
    L_outer = (L - L_active) / 2
    function reference_dt(species, E1_ref, dE_ref)
      E0_ref = E1_ref - dE_ref
      dt_ref = L_outer/E_to_v(species, E0_ref) + L_outer/E_to_v(species, E1_ref)
  
      if n_cells == 0
        L_inner = L_active / 2
        dt_ref += L_inner/E_to_v(species, E0_ref) + L_inner/E_to_v(species, E1_ref)
      else
        for i_step = 1:n_cells
          E_now_ref = E0_ref + (i_step - 1/2) * dE_ref / n_cells
          dt_ref += L_active / (n_cells * E_to_v(species, E_now_ref))
        end
      end
      return dt_ref
    end
    dt = reference_dt(species, E1_ref, dE_ref)
    dt_design = shift_time ? reference_dt(species_d, E1_d, dE_d) : dt
  else
    dt = dt_design = zero(L)
  end
  t_exit = t_enter + dt
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref(t_exit))
  else
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref)
  end
  refstate = RefState{eltype(coords.v)}(; t_enter=PreParameter(b -> b.t_ref, ()), beta_gamma_enter,
                                          t_exit=PreParameter((b, dt) -> b.t_ref + dt, (dt,)),
                                          beta_gamma_exit, L, g, ds_step,
                                          s_enter=BeamTracking.num_lower(eltype(coords.v), bunch.s),
                                          step0=bunch.step, n_steps=universalparams.n_steps)
  # As in the general universal!
  bunch.s = oftype(bunch.s, bunch.s + L)
  bunch.step = bunch.step + universalparams.n_steps
  kc = KernelChain(Val{10}(), refstate, batch_start)

  # Ramping
  if p_over_q_ref isa TimeDependentParam
    if ramp_update_each_particle
      error("ramp_update_each_particle = true not yet implemented for SaganCavity") # TODO
    end
    # As in universal! in unpack_bl.jl: the shift is exactly the change of the stored value
    p_over_q_ref_initial = bunch.p_over_q_ref
    p_over_q_ref_final = convert(typeof(bunch.p_over_q_ref), p_over_q_ref(bunch.t_ref))
    if p_over_q_ref_final != p_over_q_ref_initial
      kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (p_over_q_ref_initial, 
                                       p_over_q_ref_final-p_over_q_ref_initial, Val{!ramp_particle_energy_without_rf}())))
      setfield!(bunch, :p_over_q_ref, p_over_q_ref_final)
    end
  end

  # Entrance aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, 0, true))
        kc = @inline(aperture(kc, apertureparams, true))
      else
        kc = @inline(aperture(kc, apertureparams, true))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, 0, true))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, 0, true))
    end
  elseif isactive(apertureparams)
    kc = @inline(aperture(kc, apertureparams, true))
  end

  # Cavity tracking
  kc = @inline(sagan_cavity(tm, kc, p_over_q_ref, bunch, fullele.name, bmultipoleparams, rfparams, beamlineparams, L, ele, dE_ref, design_params))

  # Exit aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(aperture(kc, apertureparams, false))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, L, false))
      else
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, L, false))
        kc = @inline(aperture(kc, apertureparams, false))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, L, false))
    end
  elseif isactive(apertureparams)
    kc = @inline(aperture(kc, apertureparams, false))
  end

  # If the Bunch species differs from the reference species, re-reference z to the design clock
  if shift_time
    # C_LIGHT*(t_exit_design - t_exit)
    c_dt = PreParameter((b, dt, dt_design) -> C_LIGHT*((b.t_ref + dt_design) - (b.t_ref + dt)),
                        (dt, dt_design))
    kc = push(kc, make_kernel_call(BeamTracking.reference_time_shift!, (1/beta_gamma_exit, c_dt)))
  end

  # Launched by run!, which then increases bunch.t_ref by dt_design and sets bunch.p_over_q_ref: the
  # coordinates are now relative to this Beamline's reference (after the energy change dE_ref)
  p_over_q_ref_exit = p_over_q_ref isa TimeDependentParam ? bunch.p_over_q_ref : p_over_q_ref
  # As in universal! in unpack_bl.jl
  if sink isa Bunch
    @noinline run_chain!(sink, kc; kwargs...)
    sink.t_ref = sink.t_ref + dt_design
    sink.p_over_q_ref = p_over_q_ref_exit
    return nothing
  end
  return push_chain!(sink, bunch, kc, dt_design, p_over_q_ref_exit)
end

#---------------------------------------------------------------------------------------------------

# dE_ref is the reference energy change of the Bunch species, and design_params is nothing, or the 
# mass, initial reference momentum and reference energy change of the reference species if it 
# differs from the Bunch species (see universal!)
@inline function sagan_cavity(tm::SaganCavity, kc, p_over_q_ref, bunch, ele_name, bmultipoleP, rfP, beamlineP, L, ele, dE_ref, design_params)
  species = bunch.species
  mass = massof(species)
  q = chargeof(species)

  p1_over_q_ref = beamlineP.beamline.p_over_q_ref
  rf_omega = rf_omega_calc(rfP, beamlineP, ele.ref)
  E1_ref = R_to_E(species, p1_over_q_ref)
  E0_ref = E1_ref - dE_ref
  P0c = sqrt(E0_ref^2 - mass^2)
  p0q = E_to_R(species, E0_ref)
  t_phi0 = rf_phi0_calc(rfP, beamlineP.beamline.species_ref) / rf_omega

  n_cells, L_active = rf_step_calc(tm.n_cells, tm.L_active, rf_omega, L)
  L_active <= L * (1 + eps(L)) || error("Cavity $ele_name, cannot have L_active ($L_active) greater than L ($L)." )
  t_ref = 0    # bunch.t_ref ## Relative time tracking assumed for now.
  a = gyromagnetic_anomaly(species)

  #

  if L == 0
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L!,
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0, design_params)))
    else
      m_order = bmultipoleP.order
      KnL, KsL = get_integrated_strengths(bmultipoleP, L, p0q, ele.ref.field_scale)
      has_mult = (length(m_order) > 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L!, 
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), m_order, KnL.*p0q, KsL.*p0q, a, q*rfP.voltage, rf_omega, t_phi0, design_params)))
    end

  elseif L_active == 0
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L_active!,
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0, L, design_params)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, ele.ref.field_scale)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L_active!, 
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a, q*rfP.voltage, rf_omega, t_phi0, L, design_params)))
    end

  else
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], 
                a, q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L, design_params)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, ele.ref.field_scale)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a,
                q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L, design_params)))
    end
  end
end
