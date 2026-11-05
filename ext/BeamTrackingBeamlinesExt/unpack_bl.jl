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
  ele::LatElement, 
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
  # The Bunch species differing from the reference species is static for the whole track! call
  shift_time = bunch.species != ref.species_ref
  # Function barrier
  @noinline universal!(bunch, tm, ele, loweredele, p_over_q_ref, Val(ramp_particle_energy_without_rf), Val(ramp_update_each_particle), Val(shift_time), batch_start; kwargs...)
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
  ::Val{shift_time},
  batch_start;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle, shift_time}
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
  c_dt_exit = BeamTracking.num_lower(T, C_LIGHT*(t_exit_design - t_exit))
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
    kc = push(kc, make_kernel_call(BeamTracking.reference_time_shift!, (1/beta_gamma_exit, c_dt_exit)))
  end
  
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(bunch.coords, kc; kwargs...)

  # Update reference time
  bunch.t_ref = t_exit_design
  bunch.p_over_q_ref = p_over_q_ref_exit

  return nothing
end