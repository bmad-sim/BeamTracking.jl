#---------------------------------------------------------------------------------------------------
# universal! for SaganCavity tracking.
# its own universal! is needed because it does its own thing to the reference time
function universal!(
  bunch,
  tm::SaganCavity,
  fullele,
  ele,
  p_over_q_ref,
  ::Val{ramp_particle_energy_without_rf}, 
  ::Val{ramp_update_each_particle},
  ::Val{shift_time},
  batch_start;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle, shift_time}
  alignmentparams= ele.alignmentparams
  bendparams = ele.bendparams 
  bmultipoleparams = ele.bmultipoleparams 
  patchparams = ele.patchparams 
  apertureparams = ele.apertureparams 
  rfparams = ele.rfparams 
  beamlineparams = ele.beamlineparams 
  mapparams = ele.mapparams 
  fourpotentialparams = ele.fourpotentialparams 
  emultipoleparams = ele.emultipoleparams 
  emfieldparams = ele.emfieldparams 
  ibsparams= ele.ibsparams
  L = ele.universalparams.L
  coords = bunch.coords
  !shift_time || error("SaganCavity Tracking is not compatible with Bunch species differing from Beamline reference species")
  !isactive(mapparams) || error("SaganCavity Tracking through element $(fullele.name) with MapParams is undefined")
  !isactive(patchparams) || error("SaganCavity Tracking through element $(fullele.name) with PatchParams is undefined")
  !isactive(bendparams) || error("SaganCavity Tracking through element $(fullele.name) with BendParams is undefined")
  !isactive(fourpotentialparams) || error("SaganCavity Tracking through element $(fullele.name) with FourPotentialParams is undefined")
  !isactive(emultipoleparams) || error("SaganCavity Tracking through element $(fullele.name) with EMultipoleParams is undefined")
  isactive(rfparams) || error("SaganCavity Tracking through element $(fullele.name) without RFParams is undefined")

  p_over_q_ref = ele.ref.p_over_q_ref
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
  # reference time change
  if L != 0
    species = bunch.species
    p1_over_q_ref = p_over_q_ref
    rf_omega = rf_omega_calc(rfparams, beamlineparams, ele.ref)
    n_cells, L_active = rf_step_calc(tm.n_cells, tm.L_active, rf_omega, L)
    L_outer = (L - L_active) / 2
    E1_ref = R_to_E(species, p1_over_q_ref)
    dE_ref = beamlineparams.dE_ref
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
    t_exit = t_enter + dt_ref
  else
    t_exit = t_enter
  end
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref(t_exit))
  else
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref)
  end
  kc = KernelChain(Val{10}(), RefState{eltype(coords.v)}(; t_enter, beta_gamma_enter, t_exit, beta_gamma_exit, L, g, ds_step), batch_start)

  # Ramping
  if p_over_q_ref isa TimeDependentParam
    if ramp_update_each_particle
      error("ramp_update_each_particle = true not yet implemented for SaganCavity") # TODO
    end
    p_over_q_ref_initial = bunch.p_over_q_ref
    p_over_q_ref_final = p_over_q_ref(bunch.t_ref)
    if !(p_over_q_ref_initial ≈ p_over_q_ref_final)
      kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (p_over_q_ref_initial, 
                                       p_over_q_ref_final-p_over_q_ref_initial, Val{!ramp_particle_energy_without_rf}())))
      setfield!(bunch, :p_over_q_ref, p_over_q_ref_final)
    end
  end

  # Entrance aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, true))
        kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, true))
      else
        kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, true))
        kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, true))
      end
    else
      kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, true))
    end
  elseif isactive(apertureparams)
    kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, true))
  end

  # Cavity tracking
  kc = @inline(sagan_cavity(tm, kc, p_over_q_ref, bunch, fullele.name, bmultipoleparams, rfparams, beamlineparams, L, ele))

  # Exit aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, false))
        kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, false))
      else
        kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, false))
        kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, false))
      end
    else
      kc = @inline(alignment(tm, kc, p_over_q_ref, bunch, alignmentparams, bendparams, L, false))
    end
  elseif isactive(apertureparams)
    kc = @inline(aperture(tm, kc, p_over_q_ref, bunch, apertureparams, false))
  end

  # noinline necessary here for small binaries and faster execution
  @noinline launch!(coords, kc; kwargs...)
  
  bunch.t_ref = t_exit
  
  return nothing
end

#---------------------------------------------------------------------------------------------------

@inline function sagan_cavity(tm::SaganCavity, kc, p_over_q_ref, bunch, ele_name, bmultipoleP, rfP, beamlineP, L, ele)
  species = bunch.species
  mass = massof(species)
  q = chargeof(species)

  p1_over_q_ref = beamlineP.beamline.p_over_q_ref
  rf_omega = rf_omega_calc(rfP, beamlineP, ele.ref)
  E1_ref = R_to_E(species, p1_over_q_ref)
  dE_ref = beamlineP.dE_ref
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
                   t_ref, Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0)))
    else
      m_order = bmultipoleP.order
      KnL, KsL = get_integrated_strengths(bmultipoleP, L, p0q, 1)
      has_mult = (length(m_order) > 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L!, 
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), m_order, KnL.*p0q, KsL.*p0q, a, q*rfP.voltage, rf_omega, t_phi0)))
    end

  elseif L_active == 0
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L_active!,
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0, L)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, 1)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L_active!, 
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a, q*rfP.voltage, rf_omega, t_phi0, L)))
    end

  else
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], 
                a, q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, 1)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a,
                q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L)))
    end
  end
end
