# Slice k of an element (k = 1 of 1 for an element with n_slices == 1). Each slice is tracked by its
# own universal! (and launch!) like an element of length L_slice = L/n_slices starting at s0, so the
# reference energy is ramped and collective effects are recomputed at the start of each slice. L is
# still the element length (strengths are per element length), n_steps and ds_step are per slice,
# and fringe_at (F) only includes the element entrance in the first slice and the element exit in
# the last slice. fringe_at, rad_damping, rad_fluctuations, and FIRST and LAST (the first_slice and
# last_slice properties) are type parameters, and getproperty returns those, so that they are known
# to the compiler after the function barrier.
struct _SliceUniversalParams{A,F,RD,RF,FIRST,LAST}
  L::A
  n_steps::Int
  ds_step::A
  n_slices::Int
  s0::A
  L_slice::A
  # With the type parameters given (see _whole_element)
  _SliceUniversalParams{A,F,RD,RF,FIRST,LAST}(L, n_steps, ds_step, n_slices, s0, L_slice) where {A,F,RD,RF,FIRST,LAST} =
    new{A,F,RD,RF,FIRST,LAST}(L, n_steps, ds_step, n_slices, s0, L_slice)
  function _SliceUniversalParams(up::UniversalParams, k::Int)
    L = float(up.L)
    n_steps = up.n_steps
    n_slices = up.n_slices
    n_slices == 1 && return _whole_element(L, n_steps, up.fringe_at, up.rad_damping, up.rad_fluctuations)
    L_slice = L / n_slices
    ds_step = L_slice / n_steps
    s0 = (k - 1)*L_slice
    first_slice = k == 1
    last_slice = k == n_slices
    fringe_in = first_slice && (up.fringe_at == Fringe.BothEnds || up.fringe_at == Fringe.EntranceEnd)
    fringe_out = last_slice && (up.fringe_at == Fringe.BothEnds || up.fringe_at == Fringe.ExitEnd)
    fringe_at = fringe_in ? (fringe_out ? Fringe.BothEnds : Fringe.EntranceEnd) : (fringe_out ? Fringe.ExitEnd : Fringe.NoEnd)
    return new{typeof(L),fringe_at,up.rad_damping,up.rad_fluctuations,first_slice,last_slice}(L, n_steps, ds_step, n_slices, s0, L_slice)
  end
end

# The whole element (k = 1 of 1, n_slices == 1), with its type parameters compile-time constants in
# each branch. Built from the run-time values of fringe_at, rad_damping and rad_fluctuations, the
# type would be looked up in Julia's type cache at every element (~150 ns per element vs ~120 ns
# for _ConcreteUniversalParams, which it replaces, with two fewer type parameters). The most common
# case (both fringes, no radiation) is tested first
@eval function _whole_element(L::A, n_steps, fringe_at, rad_damping, rad_fluctuations) where {A}
  $(let ex = :(error("Unknown fringe_at = ", fringe_at))
      for F in (Fringe.ExitEnd, Fringe.EntranceEnd, Fringe.NoEnd, Fringe.BothEnds), RD in (true, false), RF in (true, false)
        ex = :(fringe_at == $F && rad_damping == $RD && rad_fluctuations == $RF ?
               _SliceUniversalParams{A,$F,$RD,$RF,true,true}(L, n_steps, L / n_steps, 1, zero(L), L) : $ex)
      end
      ex
    end)
end

# The first slice starts at s0 = 0 and a single slice is the whole element (L_slice = L): returned
# as such, so that an element with n_slices == 1 sees the constant 0 and L itself (not 0.0, or L/1, a
# new TPS for a TPS length)
Base.getproperty(up::_SliceUniversalParams{A,F,RD,RF,FIRST,LAST}, key::Symbol) where {A,F,RD,RF,FIRST,LAST} =
  key === :fringe_at ? F : key === :rad_damping ? RD : key === :rad_fluctuations ? RF :
  key === :first_slice ? FIRST : key === :last_slice ? LAST : key === :s0 ? (FIRST ? 0 : getfield(up, :s0)) :
  key === :L_slice ? (FIRST && LAST ? getfield(up, :L) : getfield(up, :L_slice)) : getfield(up, key)

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
  batch_start,
  pos::Base.RefValue, # position of the element in the Beamline (see universal!)
  sink;
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
  # Each slice is tracked separately (see _SliceUniversalParams)
  n_slices = (universalparams::UniversalParams).n_slices
  if n_slices > 1
    !(tm isa SaganCavity) || error("SaganCavity is not compatible with elements with n_slices > 1. Try n_slices = 1 instead")
    isnothing(patchparams) || error("Elements with n_slices > 1 are not compatible with PatchParams. Try n_slices = 1 instead")
    isnothing(mapparams) || error("Elements with n_slices > 1 are not compatible with MapParams. Try n_slices = 1 instead")
  end
  for k in 1:n_slices
    @noinline universal!(bunch, tm, ele, _SliceUniversalParams(universalparams, k), alignmentparams, bendparams, bmultipoleparams,
                         patchparams, apertureparams, rfparams, beamlineparams, mapparams, fourpotentialparams, emultipoleparams,
                         emfieldparams, ibsparams, ref, ramp_particle_energy_without_rf, ramp_update_each_particle,
                         shift_time, batch_start, pos, sink; kwargs...)
  end
end

# Step 2: Build the kernel chain and pass it to sink (see prepared_bl.jl) ----
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
  batch_start,
  pos, # Ref of the position of the slice in the Beamline (step = steps before it, s = its s), 
       # advanced past it (in a Ref: numbers passed to or returned from these dynamic calls would be
       # allocated; s has the type of bunch.t_ref, which holds the type of the lengths)
  sink;
  kwargs...
) where {ramp_particle_energy_without_rf, ramp_update_each_particle, shift_time}
  # Every element in a Beamline has the same p_over_q_ref (possibly a function of time)
  p_over_q_ref = ref.p_over_q_ref

  L = universalparams.L
  # The slice of the element tracked here (the whole element if n_slices == 1)
  s0 = universalparams.s0
  L_slice = universalparams.L_slice

  # When preparing (sink is not the Bunch tracked), what depends on the Bunch at the call is a
  # PreParameter (see prepared_bl.jl): bunch.t_ref, the reference momentum of the Bunch, and with a
  # time dependent reference energy the reference momentum at the element and everything computed
  # from it (e.g. every strength without the per-particle ramp)
  preparing = !(sink isa Bunch)
  ibs = !isnothing(ibsparams) && ibsparams.ibs_model != IBSModel.None && L_slice > 0
  ramp_per_particle = p_over_q_ref isa TimeDependentParam && ramp_update_each_particle

  # Compute information about reference coordinate system:
  t_enter = bunch.t_ref
  # Evaluate p_over_q_ref first and convert the number: TimeDependentParams of the energy are 
  # only needed for per-particle ramping
  if p_over_q_ref isa TimeDependentParam
    p_over_q_ref_enter = preparing ? PreParameter((b, p) -> p(b.t_ref), (p_over_q_ref,)) : p_over_q_ref(t_enter)
  else
    p_over_q_ref_enter = p_over_q_ref
  end
  beta_gamma_enter = R_to_beta_gamma(bunch.species, p_over_q_ref_enter)
  g = isnothing(bendparams) ? (0,0) : reverse((bendparams.g_ref .* sincos(bendparams.tilt_ref)))
  ds_step = universalparams.ds_step

  # Reference time evolution thru element assumes constant energy
  # using the energy at the start of the element.
  # bunch.t_ref is the design clock: the time of the Beamline reference species at the Beamline
  # reference energy, which time dependent parameters and RF are timed against. t_exit is the time
  # of the Bunch reference particle, which the kernels' z is measured against. The two differ only
  # if the Bunch species differs from the reference species, in which case z is re-referenced to 
  # the design clock at the end of the element (slice).
  dt = L_slice / beta_gamma_to_v(beta_gamma_enter)
  dt_design = shift_time ? L_slice / R_to_v(ref.species_ref, ref.field_scale*p_over_q_ref_enter) : dt

  if ramp_per_particle
    # In the number type of bunch.p_over_q_ref, so that the particles are ramped to exactly the
    # reference momentum the Bunch stores after the element
    if preparing
      p_over_q_ref_exit = PreParameter((b, p, dtd) -> convert(typeof(b.p_over_q_ref), p(b.t_ref + dtd)),
                                       (p_over_q_ref, dt_design))
    else
      p_over_q_ref_exit = convert(typeof(bunch.p_over_q_ref), p_over_q_ref(bunch.t_ref + dt_design))
    end
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
    t_enter = PreParameter(b -> BeamTracking.num_lower(eltype(b.coords.v), b.t_ref), ()),
    beta_gamma_enter = BeamTracking.num_lower(T, beta_gamma_enter),
    t_exit = PreParameter((b, dt) -> BeamTracking.num_lower(eltype(b.coords.v), b.t_ref + dt),
                          (dt,)),
    beta_gamma_exit = BeamTracking.num_lower(T, beta_gamma_exit),
    L = BeamTracking.num_lower(T, L_slice),
    g = BeamTracking.num_lower(T, g),
    ds_step = BeamTracking.num_lower(T, ds_step),
    s0 = BeamTracking.num_lower(T, s0),
    s_enter = BeamTracking.num_lower(T, pos[].s),
    step0 = pos[].step,
    n_steps = universalparams.n_steps,
  )
  pos[] = (step = pos[].step + universalparams.n_steps, s = oftype(pos[].s, pos[].s + L_slice))

  # Current KernelChain length is 11 because we have up to
  # 2 aperture, 2 alignment, 1 body kernel, 1 IBS kernel,
  # 2 kernels to update the particles' reference energy,
  # 2 for coordinate conversion with implicit,
  # and 1 to update the particles' reference time
  kc = KernelChain(Val{11}(), refstate, batch_start)
  
  if ibs
    kc = @inline(ibs_kick(kc, bunch, ibsparams, bendparams, L_slice, preparing))
  end

  # With a time dependent reference energy, the reference momentum of the Bunch at the element is
  # only known at the call (it is that after the previous element)
  if preparing && p_over_q_ref isa TimeDependentParam
    bunch_beta_gamma = PreParameter(b -> R_to_beta_gamma(b.species, b.p_over_q_ref), ())
  else
    bunch_beta_gamma = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)
  end
  # Reference energy as a function of time, evaluated per particle in the ramp kernels
  beta_gamma_enter_t = ramp_per_particle ? R_to_beta_gamma(bunch.species, p_over_q_ref) : nothing

  if ramp_per_particle
    kc = push(kc,make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, beta_gamma_enter_t-bunch_beta_gamma, Val{!ramp_particle_energy_without_rf}())))
    kc = push_transforms_out(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_out! not supported yet with ramp_update_each_particle = true")))
    kc = push_transforms_in(kc, make_kernel_call((i, coords, cur_s, cur_t_ref)->error("transforms_in! not supported yet with ramp_update_each_particle = true")))
  else
    # Make sure to evaluate p_over_q_ref if not ramp_update_each_particle
    p_over_q_ref = p_over_q_ref_enter
    # bunch.p_over_q_ref becomes p_over_q_ref_enter in its own number type. The particles are shifted
    # by exactly the change of that stored value, so that their coordinates are always relative to
    # it: an approximate comparison here skipped slow ramps (below 1.5e-8 per element, 3.4e-4 with a
    # Float32 reference) but stored them anyway, so they were lost
    if p_over_q_ref isa PreParameter
      # Prepared with a time dependent reference energy: whether it changes is only known at the call,
      # so the shift is always there (reference_momentum_shift! skips a zero shift)
      p_over_q_ref_new = PreParameter((b, p) -> convert(typeof(b.p_over_q_ref), p(b.t_ref)), (ref.p_over_q_ref,))
      dbeta_gamma = R_to_beta_gamma(bunch.species, p_over_q_ref_new) - bunch_beta_gamma
      kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, dbeta_gamma, Val{!ramp_particle_energy_without_rf}())))
    else
      p_over_q_ref_new = convert(typeof(bunch.p_over_q_ref), p_over_q_ref_enter)
      if p_over_q_ref_new != bunch.p_over_q_ref
        dbeta_gamma = R_to_beta_gamma(bunch.species, p_over_q_ref_new) - bunch_beta_gamma
        kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch_beta_gamma, dbeta_gamma, Val{!ramp_particle_energy_without_rf}())))
      end
    end
  end

  # Entrance aperture (only at the element entrance) and alignment (at the slice entrance)
  entrance_apertureparams = universalparams.first_slice ? apertureparams : nothing
  if !isnothing(alignmentparams)
    if !isnothing(entrance_apertureparams)
      if entrance_apertureparams.aperture_shifts_with_body
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s0, true))
        kc = @inline(aperture(kc, entrance_apertureparams, true))
      else
        kc = @inline(aperture(kc, entrance_apertureparams, true))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s0, true))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, s0, true))
    end
  elseif !isnothing(entrance_apertureparams)
    kc = @inline(aperture(kc, entrance_apertureparams, true))
  end

  # ele.ref has the tracking p_over_q_ref of the element body: evaluated at the entrance, except 
  # with the per-particle ramp
  ele = LoweredElement(universalparams, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams, rfparams, 
                       beamlineparams, mapparams, fourpotentialparams, emultipoleparams, emfieldparams, ibsparams, 
                       TrackReference(ref.species_ref, p_over_q_ref, ref.field_scale))
  kc = @inline fringe_body_fringe(tm, kc, bunch, ele)

  # Exit aperture (only at the element exit) and alignment (at the slice exit)
  exit_apertureparams = universalparams.last_slice ? apertureparams : nothing
  s_exit = universalparams.last_slice ? L : s0 + L_slice
  if !isnothing(alignmentparams)
    if !isnothing(exit_apertureparams)
      if exit_apertureparams.aperture_shifts_with_body
        kc = @inline(aperture(kc, exit_apertureparams, false))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_exit, false))
      else
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_exit, false))
        kc = @inline(aperture(kc, exit_apertureparams, false))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_exit, false))
    end
  elseif !isnothing(exit_apertureparams)
    kc = @inline(aperture(kc, exit_apertureparams, false))
  end

  # If ramping, now need to uniformly ramp all particles to same reference energy
  if ramp_per_particle
    # If TimeDependentParam, at end ramp all 
    # uniformly to p_over_q_ref(t_ref at end)
    kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (beta_gamma_enter_t, beta_gamma_exit-beta_gamma_enter_t, Val{!ramp_particle_energy_without_rf}())))
  end

  # If the Bunch species differs from the reference species, re-reference z to the design clock
  if shift_time
    # C_LIGHT*(t_exit_design - t_exit)
    c_dt = PreParameter((b, dt, dt_design) -> C_LIGHT*((b.t_ref + dt_design) - (b.t_ref + dt)),
                        (dt, dt_design))
    kc = push(kc, make_kernel_call(BeamTracking.reference_time_shift!, (1/beta_gamma_exit, c_dt)))
  end
  
  # track! launches the chain at once, prepare stores it. The type of kc is often not inferred (it
  # depends on values, e.g. whether the reference momentum changes): it is passed on alone, as a
  # PreKernelChain it would be boxed again at every call (0.6 ms more per ESR turn)
  if sink isa Bunch
    # (with the per-particle ramp these depend on bunch.t_ref before it changes)
    dt_now = process_preparameter(dt_design, sink)
    p_over_q_ref_now = process_preparameter(p_over_q_ref_exit, sink)
    @noinline run_chain!(sink, kc; kwargs...)
    sink.t_ref = sink.t_ref + dt_now
    sink.p_over_q_ref = p_over_q_ref_now
    return nothing
  end
  return push_chain!(sink, bunch, PreKernelChain(kc, dt_design, p_over_q_ref_exit))
end