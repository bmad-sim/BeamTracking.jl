struct _ConcreteUniversalParams{A}
  L::A
  n_steps::Int
  ds_step::Float64
  n_slices::Int
  fringe_at::Fringe.T
  rad_damping::Bool
  rad_fluctuations::Bool
  function _ConcreteUniversalParams(up::UniversalParams)
    L = float(up.L)
    n_steps = up.n_steps
    ds_step = L / n_steps # actual ds_step, not fake one
    n_slices = up.n_slices
    fringe_at = up.fringe_at
    rad_damping = up.rad_damping
    rad_fluctuations = up.rad_fluctuations
    return new{typeof(L)}(L, n_steps, ds_step, n_slices, fringe_at, rad_damping, rad_fluctuations)
  end
end

Base.@kwdef struct LoweredElement{A,B,C,D,E,F,G,H,I,J,K,L,M,N}
  universalparams::A
  alignmentparams::B
  bendparams::C
  bmultipoleparams::D
  patchparams::E
  apertureparams::F
  rfparams::G
  beamlineparams::H
  mapparams::I
  fourpotentialparams::J
  emultipoleparams::K
  emfieldparams::L
  ibsparams::M
  ref::N
end

# Step 1: Unpack the element ---------------------------------------------
function _track!(
  bunch::Bunch,
  ele::LineElement, 
  context::Context,
  ref::TrackReference,
  tm,
  scalar_params,
  ramp_particle_energy_without_rf,
  ramp_update_each_particle,
  rf_on,
  batch_start;
  kwargs...
)
  # Every element in a Beamline has the same p_over_q_ref (possibly a function of time)
  p_over_q_ref = ref.p_over_q_ref

  # Parameter groups in ignore_params are replaced with nothing, exactly as if the element
  # did not have them (like RFParams with rf_on = false)
  ignoreparams = ele.IgnoreParams
  ig = isnothing(ignoreparams) ? () : ignoreparams.ignore_params
  universalparams = UniversalParams in ig ? nothing : deval(ele.UniversalParams, context)
  alignmentparams = AlignmentParams in ig ? nothing : deval(ele.AlignmentParams, context)
  bendparams = BendParams in ig ? nothing : deval(ele.BendParams, context)
  bmultipoleparams = BMultipoleParams in ig ? nothing : deval(ele.BMultipoleParams, context)
  patchparams = PatchParams in ig ? nothing : deval(ele.PatchParams, context)
  apertureparams = ApertureParams in ig ? nothing : deval(ele.ApertureParams, context)
  mapparams = MapParams in ig ? nothing : deval(ele.MapParams, context)
  rfparams = (!rf_on || RFParams in ig) ? nothing : deval(ele.RFParams, context)
  beamlineparams = BeamlineParams in ig ? nothing : deval(ele.BeamlineParams, context)
  fourpotentialparams = FourPotentialParams in ig ? nothing : deval(ele.FourPotentialParams, context)
  emultipoleparams = EMultipoleParams in ig ? nothing : deval(ele.EMultipoleParams, context)
  emfieldparams = EMFieldParams in ig ? nothing : deval(ele.EMFieldParams, context)
  ibsparams = IBSParams in ig ? nothing : deval(ele.IBSParams, context)

  universalparams = isactive(universalparams) ? universalparams : nothing
  alignmentparams = isactive(alignmentparams) ? alignmentparams : nothing
  bendparams = isactive(bendparams) ? bendparams : nothing
  bmultipoleparams = isactive(bmultipoleparams) ? bmultipoleparams : nothing
  patchparams = isactive(patchparams) ? patchparams : nothing
  apertureparams = isactive(apertureparams) ? apertureparams : nothing
  mapparams = isactive(mapparams) ? mapparams : nothing
  rfparams = isactive(rfparams) ? rfparams : nothing
  beamlineparams = isactive(beamlineparams) ? beamlineparams : nothing
  fourpotentialparams = isactive(fourpotentialparams) ? fourpotentialparams : nothing
  emultipoleparams = isactive(emultipoleparams) ? emultipoleparams : nothing
  emfieldparams = isactive(emfieldparams) ? emfieldparams : nothing
  ibsparams = isactive(ibsparams) ? ibsparams : nothing

  if scalar_params
    universalparams = scalarize(universalparams)
    alignmentparams = scalarize(alignmentparams)
    bendparams = scalarize(bendparams)
    bmultipoleparams = scalarize(bmultipoleparams)
    patchparams = scalarize(patchparams)
    apertureparams = scalarize(apertureparams)
    mapparams = scalarize(mapparams)
    rfparams = scalarize(rfparams)
    beamlineparams = scalarize(beamlineparams)
    fourpotentialparams = scalarize(fourpotentialparams)
    emultipoleparams = scalarize(emultipoleparams)
    emfieldparams = scalarize(emfieldparams)
    ibsparams = scalarize(ibsparams)
    p_over_q_ref = scalarize(p_over_q_ref)
  end

  loweredele = LoweredElement(; 
    universalparams=_ConcreteUniversalParams(universalparams),
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
  )
  # Function barrier
  @noinline universal!(bunch, tm, ele, loweredele, p_over_q_ref, Val(ramp_particle_energy_without_rf), Val(ramp_update_each_particle), batch_start; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  bunch,
  tm,
  fullele,
  ele,
  p_over_q_ref,
  ::Val{ramp_particle_energy_without_rf}, 
  ::Val{ramp_update_each_particle},
  batch_start;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle}
  universalparams = ele.universalparams
  alignmentparams = ele.alignmentparams
  bendparams = ele.bendparams
  apertureparams = ele.apertureparams
  ibsparams = ele.ibsparams

  L = universalparams.L

  if universalparams.n_slices != 1
    error("Elements with n_slices > 1 not yet supported")
  end

  # Compute information about reference coordinate system:
  t_enter = bunch.t_ref
  beta_gamma_enter_t = R_to_beta_gamma(bunch.species, p_over_q_ref)
  beta_gamma_enter = p_over_q_ref isa TimeDependentParam ? beta_gamma_enter_t(t_enter) : beta_gamma_enter_t
  g = isnothing(bendparams) ? (0,0) : reverse((bendparams.g_ref .* sincos(bendparams.tilt_ref)))
  ds_step = universalparams.ds_step

  # Reference time evolution thru element assumes constant energy
  # using the energy at the start of the element.
  # bunch.t_ref is the design clock: the time of the Beamline reference species at the Beamline
  # reference energy, which time dependent parameters and RF are timed against. t_exit is the time
  # of the Bunch reference particle, which the kernels' z is measured against. The two differ only
  # if the Bunch species differs from the reference species, in which case z is re-referenced to 
  # the design clock at the end of the element.
  t_exit = bunch.t_ref + L / beta_gamma_to_v(beta_gamma_enter)
  shift_time = bunch.species != ele.ref.species_ref
  if shift_time
    p_over_q_ref_enter = p_over_q_ref isa TimeDependentParam ? p_over_q_ref(t_enter) : p_over_q_ref
    t_exit_design = bunch.t_ref + L / R_to_v(ele.ref.species_ref, p_over_q_ref_enter)
  else
    t_exit_design = t_exit
  end
  beta_gamma_exit_t = R_to_beta_gamma(bunch.species, p_over_q_ref)

  if p_over_q_ref isa TimeDependentParam
    if ramp_update_each_particle
      beta_gamma_exit = beta_gamma_exit_t(t_exit_design)
      p_over_q_ref_exit = p_over_q_ref(t_exit_design)
    else
      beta_gamma_exit = beta_gamma_enter # Don't ramp at end of element if !ramp_update_each_particle
      p_over_q_ref_exit = p_over_q_ref(t_enter)
    end
  else
    beta_gamma_exit = beta_gamma_exit_t
    p_over_q_ref_exit = p_over_q_ref
  end

  # Lower RefState values to coordinates number types
  T = eltype(bunch.coords.v)
  t_exit_design = BeamTracking.num_lower(T, t_exit_design)
  t_enter = BeamTracking.num_lower(T, t_enter)
  beta_gamma_enter = BeamTracking.num_lower(T, beta_gamma_enter)
  t_exit  = BeamTracking.num_lower(T, t_exit)
  beta_gamma_exit  = BeamTracking.num_lower(T, beta_gamma_exit)
  L = BeamTracking.num_lower(T, L)
  g = BeamTracking.num_lower(T, g)
  ds_step = BeamTracking.num_lower(T, ds_step)
  p_over_q_ref = BeamTracking.num_lower(T, p_over_q_ref)

  # Current KernelChain length is 11 because we have up to
  # 2 aperture, 2 alignment, 1 body kernel, 1 IBS kernel,
  # 2 kernels to update the particles' reference energy,
  # 2 for coordinate conversion with implicit,
  # and 1 to update the particles' reference time
  kc = KernelChain(Val{11}(), RefState{T}(; t_enter, beta_gamma_enter, t_exit, beta_gamma_exit, L, g, ds_step), batch_start)
  
  if !isnothing(ibsparams) && ibsparams.ibs_model != IBSModel.None && L > 0
    kc = @inline(ibs_kick(kc, bunch, ibsparams, bendparams, L))
  end

  ramp_per_particle = p_over_q_ref isa TimeDependentParam && ramp_update_each_particle
  bunch_beta_gamma = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)

  if ramp_per_particle
    kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, beta_gamma_enter_t-bunch_beta_gamma, Val{!ramp_particle_energy_without_rf}())))
    kc = push_transforms_out(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_out! not supported yet with ramp_update_each_particle = true")))
    kc = push_transforms_in(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_in! not supported yet with ramp_update_each_particle = true")))
  else
    # Make sure to evaluate p_over_q_ref if not ramp_update_each_particle
    p_over_q_ref = p_over_q_ref isa TimeDependentParam ? p_over_q_ref(t_enter) : p_over_q_ref
    if !(BeamTracking.num_lower(typeof(bunch.p_over_q_ref), beta_gamma_enter) ≈ BeamTracking.num_lower(typeof(bunch.p_over_q_ref), bunch_beta_gamma))
        kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, beta_gamma_enter - bunch_beta_gamma, Val{!ramp_particle_energy_without_rf}())))
    end
  end

  # Entrance aperture and alignment
  if !isnothing(alignmentparams)
    if !isnothing(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, true))
        kc = @inline(aperture(kc, apertureparams, true))
      else
        kc = @inline(aperture(kc, apertureparams, true))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, true))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, true))
    end
  elseif !isnothing(apertureparams)
    kc = @inline(aperture(kc, apertureparams, true))
  end

  kc = @inline fringe_body_fringe(tm, kc, p_over_q_ref, bunch, ele)

  # Exit aperture and alignment
  if !isnothing(alignmentparams)
    if !isnothing(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = @inline(aperture(kc, apertureparams, false))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, false))
      else
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, false))
        kc = @inline(aperture(kc, apertureparams, false))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, false))
    end
  elseif !isnothing(apertureparams)
    kc = @inline(aperture(kc, apertureparams, false))
  end

  # If ramping, now need to uniformly ramp all particles to same reference energy
  if ramp_per_particle
    # If TimeDependentParam, at end ramp all 
    # uniformly to p_over_q_ref(t_ref at end)
    kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (beta_gamma_enter_t, beta_gamma_exit-beta_gamma_enter_t, Val{!ramp_particle_energy_without_rf}())))
  end

  # If the Bunch species differs from the reference species, re-reference z to the design clock
  if shift_time
    kc = push(kc, make_kernel_call(BeamTracking.reference_time_shift!, (1/beta_gamma_exit, C_LIGHT*(t_exit_design - t_exit))))
  end
  
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(bunch.coords, kc; kwargs...)

  # Update reference time
  bunch.t_ref = t_exit_design
  bunch.p_over_q_ref = p_over_q_ref_exit

  return nothing
end

#---------------------------------------------------------------------------------------------------
# universal! for SaganCavity tracking.

function universal!(coords, tm::SaganCavity, ele, ramp_particle_energy_without_rf, ramp_update_each_particle, batch_start, bunch, L,
  p_over_q_ref, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams,
  rfparams, beamlineparams, mapparams, fourpotentialparams, emultipoleparams, em_field_params; kwargs...)

  !isactive(mapparams) || error("SaganCavity Tracking through element $(ele.name) with MapParams is undefined")
  !isactive(patchparams) || error("SaganCavity Tracking through element $(ele.name) with PatchParams is undefined")
  !isactive(bendparams) || error("SaganCavity Tracking through element $(ele.name) with BendParams is undefined")
  !isactive(fourpotentialparams) || error("SaganCavity Tracking through element $(ele.name) with FourPotentialParams is undefined")
  !isactive(emultipoleparams) || error("SaganCavity Tracking through element $(ele.name) with EMultipoleParams is undefined")
  isactive(rfparams) || error("SaganCavity Tracking through element $(ele.name) without RFParams is undefined")

  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)

  # Compute information about reference coordinate system:
  t_enter = bunch.t_ref
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref(t_enter))
  else
    beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref)
  end
  g = isnothing(bendparams) ? (0,0) : reverse((bendparams.g_ref .* sincos(bendparams.tilt_ref)))
  ds_step = (L == 0 || isactive(patchparams)) ? L : BeamTracking.find_steps(tm, L)[2]
  # reference time change
  if L != 0
    species = bunch.species
    p1_over_q_ref = p_over_q_ref
    rf_omega = rf_omega_calc(rfparams, beamlineparams)
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
  kc = @inline(sagan_cavity(tm, kc, p_over_q_ref, bunch, ele.name, bmultipoleparams, rfparams, beamlineparams, L))

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

# === Drift === #
@inline drift(tm, kc, p_over_q_ref, bunch, L) = error("Undefined for tracking method $tm")

# == Implicit === #
@inline implicit(tm, kc, p_over_q_ref, bunch, fourpotentialparams, bendparams, L) = error("Undefined for tracking method $tm")

# === Straight Elements === #
# "Pure" means only ONE SINGLE multipole.
# When "pure" is not present, it means that at least one HIGHER ORDER
# multipole exists.
@inline thin_pure_rf(tm, kc, p_over_q_ref, bunch, rfparams)                          = error("Undefined for tracking method $tm")
@inline thin_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bm0)                        = error("Undefined for tracking method $tm")
@inline thin_bsolenoid(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                = error("Undefined for tracking method $tm")
@inline thin_pure_bdipole(tm, kc, p_over_q_ref, bunch, bm1)                          = error("Undefined for tracking method $tm")
@inline thin_bdipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                  = error("Undefined for tracking method $tm")
@inline thin_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bm2)                      = error("Undefined for tracking method $tm")
@inline thin_bquadrupole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)              = error("Undefined for tracking method $tm")
@inline thin_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bmk)                       = error("Undefined for tracking method $tm")
@inline thin_bmultipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)               = error("Undefined for tracking method $tm")
@inline thin_bmultipole_rf(tm, kc, p_over_q_ref, bunch, bmultipoleparams, rfparams)  = error("Undefined for tracking method $tm")
@inline thin_pure_edipole(tm, kc, p_over_q_ref, bunch, em1)                          = error("Undefined for tracking method $tm")

@inline thick_pure_rf(tm, kc, p_over_q_ref, bunch, rfparams, beamlineparams, L)                         = error("Undefined for tracking method $tm")
@inline thick_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bm0, L)                                       = error("Undefined for tracking method $tm")
@inline thick_bsolenoid(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                               = error("Undefined for tracking method $tm")
@inline thick_pure_bdipole(tm, kc, p_over_q_ref, bunch, bm1, L)                                         = error("Undefined for tracking method $tm")
@inline thick_bdipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                                 = error("Undefined for tracking method $tm")
@inline thick_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bm2, L)                                     = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                             = error("Undefined for tracking method $tm")
@inline thick_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bmk, L)                                      = error("Undefined for tracking method $tm")
@inline thick_bmultipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                              = error("Undefined for tracking method $tm")
@inline thick_bmultipole_rf(tm, kc, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L) = error("Undefined for tracking method $tm")
@inline thick_pure_edipole(tm, kc, p_over_q_ref, bunch, em1, L)                                         = error("Undefined for tracking method $tm")

# === Elements with curving coordinate system "bend" === #
# "Bend" means ONLY a coordinate system curvature through the element.
# It does NOT IMPLY ANY DIPOLE FIELD! Bend specifies if the integration 
# path is curving but does not IMPACT THE PHYSICS inside.
# "Pure" means only ONE SINGLE MULTIPOLE
# When "pure" is not present, it means that at least one higher order 
# multipole exists.

# SciBmad will probably not support thin bends ever but I leave them here for now
@inline thin_bend_no_field(tm, kc, p_over_q_ref, bunch, bendparams)                      = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bm0)           = error("Undefined for tracking method $tm")
@inline thin_bend_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)   = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bm1)             = error("Undefined for tracking method $tm")
@inline thin_bend_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)     = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bm2)         = error("Undefined for tracking method $tm")
@inline thin_bend_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams) = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmk)          = error("Undefined for tracking method $tm")
@inline thin_bend_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)  = error("Undefined for tracking method $tm")

@inline thick_bend_no_field(tm, kc, p_over_q_ref, bunch, bendparams, L)                      = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bm0, L)           = error("Undefined for tracking method $tm")
@inline thick_bend_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)   = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bm1, L)             = error("Undefined for tracking method $tm")
@inline thick_bend_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)     = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bm2, L)         = error("Undefined for tracking method $tm")
@inline thick_bend_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L) = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmk, L)          = error("Undefined for tracking method $tm")
@inline thick_bend_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)  = error("Undefined for tracking method $tm")


# === Elements thin vs thick check === #
@inline pure_rf(tm, kc, p_over_q_ref, bunch, rfparams, beamlineparams, L)                          = L == 0 ? thin_pure_rf(tm, kc, p_over_q_ref, bunch, rfparams, beamlineparams)                         : thick_pure_rf(tm, kc, p_over_q_ref, bunch, rfparams, beamlineparams, L)
@inline pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bm0, L)                                   = L == 0 ? thin_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bm0)                                  : thick_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bm0, L)      
@inline bsolenoid(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                           = L == 0 ? thin_bsolenoid(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                          : thick_bsolenoid(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)       
@inline pure_bdipole(tm, kc, p_over_q_ref, bunch, bm1, L)                                     = L == 0 ? thin_pure_bdipole(tm, kc, p_over_q_ref, bunch, bm1)                                    : thick_pure_bdipole(tm, kc, p_over_q_ref, bunch, bm1, L)          
@inline bdipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                             = L == 0 ? thin_bdipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                            : thick_bdipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)             
@inline pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bm2, L)                                 = L == 0 ? thin_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bm2)                                : thick_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bm2, L)        
@inline bquadrupole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                         = L == 0 ? thin_bquadrupole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                        : thick_bquadrupole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)           
@inline pure_bmultipole(tm, kc, p_over_q_ref, bunch, bmk, L)                                  = L == 0 ? thin_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bmk)                                 : thick_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bmk, L)                   
@inline bmultipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)                          = L == 0 ? thin_bmultipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams)                         : thick_bmultipole(tm, kc, p_over_q_ref, bunch, bmultipoleparams, L)
@inline bmultipole_rf(tm, kc, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L)  = L == 0 ? thin_bmultipole_rf(tm, kc, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams) : thick_bmultipole_rf(tm, kc, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L)        
@inline bend_no_field(tm, kc, p_over_q_ref, bunch, bendparams, L)                             = L == 0 ? thin_bend_no_field(tm, kc, p_over_q_ref, bunch, bendparams)                            : thick_bend_no_field(tm, kc, p_over_q_ref, bunch, bendparams, L)
@inline bend_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bm0, L)                  = L == 0 ? thin_bend_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bm0)                 : thick_bend_pure_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bm0, L)      
@inline bend_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)          = L == 0 ? thin_bend_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)         : thick_bend_bsolenoid(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)         
@inline bend_pure_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bm1, L)                    = L == 0 ? thin_bend_pure_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bm1)                   : thick_bend_pure_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bm1, L)          
@inline bend_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)            = L == 0 ? thin_bend_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)           : thick_bend_bdipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)             
@inline bend_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bm2, L)                = L == 0 ? thin_bend_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bm2)               : thick_bend_pure_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bm2, L)        
@inline bend_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)        = L == 0 ? thin_bend_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)       : thick_bend_bquadrupole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)           
@inline bend_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmk, L)                 = L == 0 ? thin_bend_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmk)                : thick_bend_pure_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmk, L)                   
@inline bend_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)         = L == 0 ? thin_bend_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams)        : thick_bend_bmultipole(tm, kc, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)      
@inline pure_edipole(tm, kc, p_over_q_ref, bunch, em1, L)                                     = L == 0 ? thin_pure_edipole(tm, kc, p_over_q_ref, bunch, em1)                                    : thick_pure_edipole(tm, kc, p_over_q_ref, bunch, em1, L)                          