# fringe_at, rad_damping and rad_fluctuations are also type parameters, and getproperty returns 
# those, so that they are known to the compiler after the function barrier
struct _ConcreteUniversalParams{A,F,RD,RF}
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
    return new{typeof(L),fringe_at,rad_damping,rad_fluctuations}(L, n_steps, ds_step, n_slices, fringe_at, rad_damping, rad_fluctuations)
  end
end

Base.getproperty(up::_ConcreteUniversalParams{A,F,RD,RF}, key::Symbol) where {A,F,RD,RF} = 
  key === :fringe_at ? F : key === :rad_damping ? RD : key === :rad_fluctuations ? RF : getfield(up, key)

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
  scalar_params,
  ramp_particle_energy_without_rf::Val,
  ramp_update_each_particle::Val,
  shift_time::Val,
  rf_on,
  batch_start;
  kwargs...
)
  # Get the parameter groups in one pass over the element's ParamDict, then its parent's (if it 
  # has InheritParams), etc. A child's parameter group takes precedence over its parent's. This 
  # is the same as ele.<parameter group> for each group, but without hashing a type for each lookup
  universalparams = nothing
  alignmentparams = nothing
  bendparams = nothing
  bmultipoleparams = nothing
  patchparams = nothing
  apertureparams = nothing
  mapparams = nothing
  rfparams = nothing
  beamlineparams = nothing
  fourpotentialparams = nothing
  emultipoleparams = nothing
  emfieldparams = nothing
  ibsparams = nothing
  ignoreparams = nothing
  pdict = getfield(ele, :pdict)
  while true
    # As in getproperty, a protected parameter group errors unless a child already had it
    protect = get(pdict, Beamlines.ProtectParams, nothing)
    if !isnothing(protect)
      for (G, p) in (UniversalParams=>universalparams, AlignmentParams=>alignmentparams, BendParams=>bendparams, 
                     BMultipoleParams=>bmultipoleparams, PatchParams=>patchparams, ApertureParams=>apertureparams, 
                     MapParams=>mapparams, RFParams=>rfparams, BeamlineParams=>beamlineparams, 
                     FourPotentialParams=>fourpotentialparams, EMultipoleParams=>emultipoleparams, 
                     EMFieldParams=>emfieldparams, IBSParams=>ibsparams, IgnoreParams=>ignoreparams)
        isnothing(p) && nameof(G) in (protect::Beamlines.ProtectParams).protected_properties && error("Cannot get $G: parameter group is protected by ProtectParams. This can be unsafely-overridden using `unsafe_getparams`")
      end
    end
    for (T, p) in pdict
      if T === UniversalParams
        universalparams = isnothing(universalparams) ? p : universalparams
      elseif T === AlignmentParams
        alignmentparams = isnothing(alignmentparams) ? p : alignmentparams
      elseif T === BendParams
        bendparams = isnothing(bendparams) ? p : bendparams
      elseif T === BMultipoleParams
        bmultipoleparams = isnothing(bmultipoleparams) ? p : bmultipoleparams
      elseif T === PatchParams
        patchparams = isnothing(patchparams) ? p : patchparams
      elseif T === ApertureParams
        apertureparams = isnothing(apertureparams) ? p : apertureparams
      elseif T === MapParams
        mapparams = isnothing(mapparams) ? p : mapparams
      elseif T === RFParams
        rfparams = isnothing(rfparams) ? p : rfparams
      elseif T === BeamlineParams
        beamlineparams = isnothing(beamlineparams) ? p : beamlineparams
      elseif T === FourPotentialParams
        fourpotentialparams = isnothing(fourpotentialparams) ? p : fourpotentialparams
      elseif T === EMultipoleParams
        emultipoleparams = isnothing(emultipoleparams) ? p : emultipoleparams
      elseif T === EMFieldParams
        emfieldparams = isnothing(emfieldparams) ? p : emfieldparams
      elseif T === IBSParams
        ibsparams = isnothing(ibsparams) ? p : ibsparams
      elseif T === IgnoreParams
        ignoreparams = isnothing(ignoreparams) ? p : ignoreparams
      end
    end
    inherit = get(pdict, Beamlines.InheritParams, nothing)
    isnothing(inherit) && break
    pdict = getfield((inherit::Beamlines.InheritParams).parent::LatElement, :pdict)
  end

  # Parameter groups in ignore_params are replaced with nothing, exactly as if the element
  # did not have them (like RFParams with rf_on = false)
  ig = isnothing(ignoreparams) ? () : (ignoreparams::IgnoreParams).ignore_params
  universalparams = (UniversalParams in ig || isnothing(universalparams)) ? nothing : deval(universalparams, context)
  alignmentparams = (AlignmentParams in ig || isnothing(alignmentparams)) ? nothing : deval(alignmentparams, context)
  bendparams = (BendParams in ig || isnothing(bendparams)) ? nothing : deval(bendparams, context)
  bmultipoleparams = (BMultipoleParams in ig || isnothing(bmultipoleparams)) ? nothing : deval(bmultipoleparams, context)
  patchparams = (PatchParams in ig || isnothing(patchparams)) ? nothing : deval(patchparams, context)
  apertureparams = (ApertureParams in ig || isnothing(apertureparams)) ? nothing : deval(apertureparams, context)
  mapparams = (MapParams in ig || isnothing(mapparams)) ? nothing : deval(mapparams, context)
  rfparams = (!rf_on || RFParams in ig || isnothing(rfparams)) ? nothing : deval(rfparams, context)
  beamlineparams = (BeamlineParams in ig || isnothing(beamlineparams)) ? nothing : deval(beamlineparams, context)
  fourpotentialparams = (FourPotentialParams in ig || isnothing(fourpotentialparams)) ? nothing : deval(fourpotentialparams, context)
  emultipoleparams = (EMultipoleParams in ig || isnothing(emultipoleparams)) ? nothing : deval(emultipoleparams, context)
  emfieldparams = (EMFieldParams in ig || isnothing(emfieldparams)) ? nothing : deval(emfieldparams, context)
  ibsparams = (IBSParams in ig || isnothing(ibsparams)) ? nothing : deval(ibsparams, context)
  # Same as ele.tracking_method, without another lookup through the parent element
  tm = isnothing(universalparams) ? nothing : universalparams.tracking_method
  !isnothing(tm) || error("LatElement $(ele.name) has no tracking_method set")

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
    ref = TrackReference(ref.species_ref, scalarize(ref.p_over_q_ref), ref.field_scale)
  end

  # Function barrier. The LoweredElement is built after the barrier, where all types are known
  @noinline universal!(bunch, tm, ele, _ConcreteUniversalParams(universalparams), alignmentparams, bendparams, bmultipoleparams, 
                       patchparams, apertureparams, rfparams, beamlineparams, mapparams, fourpotentialparams, emultipoleparams, 
                       emfieldparams, ibsparams, ref, ramp_particle_energy_without_rf, ramp_update_each_particle, 
                       shift_time, batch_start; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  bunch,
  tm,
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
  batch_start;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle, shift_time}
  ele = LoweredElement(universalparams, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams, rfparams, 
                       beamlineparams, mapparams, fourpotentialparams, emultipoleparams, emfieldparams, ibsparams, ref)
  # Every element in a Beamline has the same p_over_q_ref (possibly a function of time)
  p_over_q_ref = ref.p_over_q_ref

  L = universalparams.L

  if universalparams.n_slices != 1
    error("Elements with n_slices > 1 not yet supported")
  end

  # Compute information about reference coordinate system:
  t_enter = bunch.t_ref
  # Evaluate p_over_q_ref first and convert the number: TimeDependentParams of the energy are 
  # only needed for per-particle ramping
  p_over_q_ref_enter = p_over_q_ref isa TimeDependentParam ? p_over_q_ref(t_enter) : p_over_q_ref
  beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref_enter)
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
    t_exit_design = bunch.t_ref + L / R_to_v(ele.ref.species_ref, ele.ref.field_scale*p_over_q_ref_enter)
  else
    t_exit_design = t_exit
  end

  if p_over_q_ref isa TimeDependentParam && ramp_update_each_particle
    p_over_q_ref_exit = p_over_q_ref(t_exit_design)
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref_exit)
  else
    # Don't ramp at end of element if !ramp_update_each_particle
    p_over_q_ref_exit = p_over_q_ref_enter
    beta_gamma_exit = beta_gamma_enter
  end

  # Lower RefState values to coordinates number types. Only the RefState gets the lowered values:
  # everything else below is computed at full precision (e.g. differences of nearly equal times
  # and energies, and bunch.t_ref), and push lowers the kernel arguments itself
  T = eltype(bunch.coords.v)
  refstate = RefState{T}(;
    t_enter = BeamTracking.num_lower(T, t_enter),
    beta_gamma_enter = BeamTracking.num_lower(T, beta_gamma_enter),
    t_exit = BeamTracking.num_lower(T, t_exit),
    beta_gamma_exit = BeamTracking.num_lower(T, beta_gamma_exit),
    L = BeamTracking.num_lower(T, L),
    g = BeamTracking.num_lower(T, g),
    ds_step = BeamTracking.num_lower(T, ds_step),
  )

  # Current KernelChain length is 11 because we have up to
  # 2 aperture, 2 alignment, 1 body kernel, 1 IBS kernel,
  # 2 kernels to update the particles' reference energy,
  # 2 for coordinate conversion with implicit,
  # and 1 to update the particles' reference time
  kc = KernelChain(Val{11}(), refstate, batch_start)
  
  if !isnothing(ibsparams) && ibsparams.ibs_model != IBSModel.None && L > 0
    kc = @inline(ibs_kick(kc, bunch, ibsparams, bendparams, L))
  end

  ramp_per_particle = p_over_q_ref isa TimeDependentParam && ramp_update_each_particle
  bunch_beta_gamma = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)
  # Reference energy as a function of time, evaluated per particle in the ramp kernels
  beta_gamma_enter_t = ramp_per_particle ? R_to_beta_gamma(bunch.species, p_over_q_ref) : nothing

  if ramp_per_particle
    kc = push(kc,make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, beta_gamma_enter_t-bunch_beta_gamma, Val{!ramp_particle_energy_without_rf}())))
    kc = push_transforms_out(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_out! not supported yet with ramp_update_each_particle = true")))
    kc = push_transforms_in(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_in! not supported yet with ramp_update_each_particle = true")))
  else
    # Make sure to evaluate p_over_q_ref if not ramp_update_each_particle
    p_over_q_ref = p_over_q_ref_enter
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