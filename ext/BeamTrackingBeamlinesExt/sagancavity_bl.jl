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
  # The slice [s0, s_end] of the element tracked here (the whole element if n_slices == 1)
  first_slice = universalparams.first_slice
  last_slice = universalparams.last_slice
  s0 = universalparams.s0
  L_slice = universalparams.L_slice
  sliced = !(first_slice && last_slice)
  !(sliced && p_over_q_ref isa TimeDependentParam) || error("SaganCavity with n_slices > 1 is not compatible with a time dependent reference energy. Try n_slices = 1 instead")
  # s_end must be the s0 of the next slice exactly, (k-1)*L_slice for slice k (s0 + L_slice can 
  # differ by rounding), so that every kick and fringe is in exactly one slice (see sagan_window)
  if last_slice
    s_end = L
  else
    k = count(j -> j*L_slice <= s0, 1:universalparams.n_slices-1) + 1
    s_end = k*L_slice
  end
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

  # Fringes of the element: the slice fringe_at only has those at the element ends, but these are at 
  # the ends of the active region, which can be inside any slice (see sagan_window)
  element_fringe_at = fullele.fringe_at
  fringe_on_in = element_fringe_at == Fringe.BothEnds || element_fringe_at == Fringe.EntranceEnd
  fringe_on_out = element_fringe_at == Fringe.BothEnds || element_fringe_at == Fringe.ExitEnd

  # reference time change, and reference energy at the slice exit
  if L != 0
    rf_omega = rf_omega_calc(rfparams, beamlineparams, ele.ref)
    n_cells, L_active = rf_step_calc(tm.n_cells, tm.L_active, rf_omega, L)
    window = sagan_window(n_cells, L_active, L, s0, s_end, first_slice, last_slice, fringe_on_in, fringe_on_out)
    dt, E_enter, E_exit = sagan_window_ref(species, E1_ref, dE_ref, n_cells, L_active, L, window)
    dt_design = shift_time ? sagan_window_ref(species_d, E1_d, dE_d, n_cells, L_active, L, window)[1] : dt
  else
    window = nothing
    E_enter = E1_ref - dE_ref
    E_exit = E1_ref
    dt = dt_design = zero(L)
  end
  t_exit = t_enter + dt
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref(t_exit))
  elseif last_slice
    beta_gamma_exit = R_to_beta_gamma(bunch.species, p_over_q_ref)
  else
    # Inside the cavity: the reference energy after the kicks of this and the earlier slices
    beta_gamma_exit = R_to_beta_gamma(bunch.species, E_to_R(species, E_exit))
  end
  if sliced
    beta_gamma_enter = R_to_beta_gamma(bunch.species, E_to_R(species, E_enter))
  end
  refstate = RefState{eltype(coords.v)}(; t_enter=PreParameter(b -> b.t_ref, ()), beta_gamma_enter,
                                          t_exit=PreParameter((b, dt) -> b.t_ref + dt, (dt,)),
                                          beta_gamma_exit, L = L_slice, g, ds_step, s0)
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

  # Entrance aperture (only at the element entrance) and alignment (at the slice entrance)
  entrance_apertureparams = first_slice && isactive(apertureparams) ? apertureparams : nothing
  if isactive(alignmentparams)
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

  # Cavity tracking
  kc = @inline(sagan_cavity(tm, kc, p_over_q_ref, bunch, fullele.name, bmultipoleparams, rfparams, beamlineparams, L, ele, dE_ref, design_params, window))

  # Exit aperture (only at the element exit) and alignment (at the slice exit)
  exit_apertureparams = last_slice && isactive(apertureparams) ? apertureparams : nothing
  if isactive(alignmentparams)
    if !isnothing(exit_apertureparams)
      if exit_apertureparams.aperture_shifts_with_body
        kc = @inline(aperture(kc, exit_apertureparams, false))
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_end, false))
      else
        kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_end, false))
        kc = @inline(aperture(kc, exit_apertureparams, false))
      end
    else
      kc = @inline(alignment(kc, alignmentparams, bendparams, L, s_end, false))
    end
  elseif !isnothing(exit_apertureparams)
    kc = @inline(aperture(kc, exit_apertureparams, false))
  end

  # If the Bunch species differs from the reference species, re-reference z to the design clock
  if shift_time
    # C_LIGHT*(t_exit_design - t_exit)
    c_dt = PreParameter((b, dt, dt_design) -> C_LIGHT*((b.t_ref + dt_design) - (b.t_ref + dt)),
                        (dt, dt_design))
    kc = push(kc, make_kernel_call(BeamTracking.reference_time_shift!, (1/beta_gamma_exit, c_dt)))
  end

  # Launched by run!, which then increases bunch.t_ref by dt_design and sets bunch.p_over_q_ref: the
  # coordinates are now relative to this Beamline's reference (after the energy change dE_ref), or 
  # inside the cavity to the reference energy at the slice exit, so the next slice does not shift 
  # them again
  if p_over_q_ref isa TimeDependentParam
    p_over_q_ref_exit = bunch.p_over_q_ref
  elseif last_slice
    p_over_q_ref_exit = p_over_q_ref
  else
    p_over_q_ref_exit = convert(typeof(bunch.p_over_q_ref), E_to_R(species, E_exit))
  end
  # As in universal! in unpack_bl.jl
  if sink isa Bunch
    @noinline run_chain!(sink, kc; kwargs...)
    sink.t_ref = sink.t_ref + dt_design
    sink.p_over_q_ref = p_over_q_ref_exit
    return nothing
  end
  return push_chain!(sink, bunch, PreKernelChain(kc, dt_design, p_over_q_ref_exit))
end

#---------------------------------------------------------------------------------------------------

# dE_ref is the reference energy change of the Bunch species, and design_params is nothing, or the 
# mass, initial reference momentum and reference energy change of the reference species if it 
# differs from the Bunch species (see universal!)
@inline function sagan_cavity(tm::SaganCavity, kc, p_over_q_ref, bunch, ele_name, bmultipoleP, rfP, beamlineP, L, ele, dE_ref, design_params, window)
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
                   t_ref, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], a, q*rfP.voltage, rf_omega, t_phi0, L, window, design_params)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, ele.ref.field_scale)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_zero_L_active!, 
                   (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), mass, q, P0c, dE_ref,
                   t_ref, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a, q*rfP.voltage, rf_omega, t_phi0, L, window, design_params)))
    end

  else
    if isnothing(bmultipoleP)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{false}(), Val{false}(), SA{Int32}[], SA{Int32}[], SA{Int32}[], 
                a, q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L, window, design_params)))
    else
      m_order = bmultipoleP.order
      Kn, Ks = get_strengths(bmultipoleP, L, p0q, ele.ref.field_scale)
      has_mult = (length(m_order) > 0)
      has_sol = (has_mult && m_order[1] == 0)
      return push(kc, make_kernel_call(BeamTracking.sagan_cavity_thick!, 
                (Val{ele.universalparams.rad_damping}(), Val{ele.universalparams.rad_fluctuations}(), Val{rfP.traveling_wave}(), 
                mass, q, P0c, dE_ref, t_ref, n_cells, Val{has_mult}(), Val{has_sol}(), m_order, Kn.*p0q, Ks.*p0q, a,
                q*rfP.voltage/L_active, rf_omega, t_phi0, L_active, L, window, design_params)))
    end
  end
end

#---------------------------------------------------------------------------------------------------
"""
    sagan_window(n_cells, L_active, L, s0, s_end, first_slice, last_slice, fringe_on_in, fringe_on_out)

The parts of a SaganCavity of length `L` that a slice from `s0` to `s_end` tracks (the window 
argument of the `sagan_cavity_thick!` and `sagan_cavity_zero_L_active!` kernels; see there).

Slicing never changes the cavity: its kicks (at fixed positions set by `n_cells` and `L_active`), 
fringes and drifts are tracked by the slice they are in, and a drift that crosses a cut is split (its 
half kicks stay at its ends, with the strength of the whole drift). A 
kick or fringe at a position `s` belongs to the slice with `s0 <= s < s_end` (a kick on a cut stays 
whole, in the slice that starts there); anything before 0 or after `L` by rounding belongs to the 
first or last slice. `s_end` must be exactly the `s0` of the next slice. With `s0 = 0`, `s_end = L` 
(one slice), the window is the whole cavity.
"""
function sagan_window(n_cells, L_active, L, s0, s_end, first_slice, last_slice, fringe_on_in, fringe_on_out)
  in_window(s) = (first_slice || s >= s0) && (last_slice || s < s_end)
  # Whether the part of the drift from a to b in the window is tracked, whether all of it is, its 
  # length (only used if not all of it is), and whether the half kicks at its start and end (at a and 
  # b, with the strength of the whole drift) are in the window
  function part(a, b)
    full = (first_slice || a >= s0) && (last_slice || b <= s_end)
    lo = first_slice ? a : max(a, s0)
    hi = last_slice ? b : min(b, s_end)
    len = full ? zero(hi - lo) : max(hi - lo, zero(hi - lo))
    kick_in = full || ((first_slice || a >= s0) && (last_slice || a < s_end))
    kick_out = full || ((first_slice || b > s0) && (last_slice || b <= s_end))
    return (full || len > 0, full, len, kick_in, kick_out)
  end
  no_part = (false, false, zero(L - s0), false, false)

  L_outer = (L - L_active)/2
  n_kicks = (L_active == 0 || n_cells == 0) ? 1 : n_cells + 1
  kick_s(j) = L_active == 0 ? L/2 : (n_cells == 0 ? L_outer + L_active/2 : L_outer + j*(L_active/n_cells))
  s_active_in = L_active == 0 ? L/2 : L_outer
  s_active_out = L_active == 0 ? L/2 : (n_cells == 0 ? L_outer + L_active : kick_s(n_cells))
  # Inside drift i: before kick i and after kick i-1
  function inside(i)
    if L_active == 0
      return no_part
    elseif n_cells == 0
      return i == 0 ? part(s_active_in, kick_s(0)) : part(kick_s(0), s_active_out)
    else
      return 1 <= i <= n_cells ? part(kick_s(i-1), kick_s(i)) : no_part
    end
  end

  # Kicks ka:kb (0-based) are in the window, kicks 0:ka-1 before it
  ka = first_slice ? 0 : count(j -> kick_s(j) < s0, 0:n_kicks-1)
  kb = last_slice ? n_kicks-1 : count(j -> kick_s(j) < s_end, 0:n_kicks-1) - 1

  outer_in, outer_in_full, outer_in_len, outer_in_kick_in, outer_in_kick_out = part(zero(L), s_active_in)
  pre, pre_full, pre_len, pre_kick_in, pre_kick_out = inside(ka)
  post, post_full, post_len, post_kick_in, post_kick_out = kb >= ka ? inside(kb+1) : no_part
  outer_out, outer_out_full, outer_out_len, outer_out_kick_in, outer_out_kick_out = part(s_active_out, L)
  fringe_in = L_active != 0 && fringe_on_in && in_window(s_active_in)
  fringe_out = L_active != 0 && fringe_on_out && in_window(s_active_out)
  return (; outer_in, outer_in_full, outer_in_len, outer_in_kick_in, outer_in_kick_out, fringe_in, 
            pre, pre_full, pre_len, pre_kick_in, pre_kick_out, ka, kb, 
            post, post_full, post_len, post_kick_in, post_kick_out, fringe_out, 
            outer_out, outer_out_full, outer_out_len, outer_out_kick_in, outer_out_kick_out)
end

"""
    sagan_window_ref(species, E1_ref, dE_ref, n_cells, L_active, L, window) -> dt_ref, E_enter, E_exit

The reference time of flight of `species` through the `window` of a SaganCavity (see `sagan_window`) 
whose reference energy goes from `E1_ref - dE_ref` to `E1_ref`, and the reference energies at the 
start and end of the window. For the whole cavity, `dt_ref` is the same (bitwise) as without slices.
"""
function sagan_window_ref(species, E1_ref, dE_ref, n_cells, L_active, L, w)
  E0_ref = E1_ref - dE_ref
  L_outer = (L - L_active)/2
  n_kicks = (L_active == 0 || n_cells == 0) ? 1 : n_cells + 1
  # Reference energy after kicks 0:j-1, i.e. along inside drift j
  E_after(j) = j <= 0 ? E0_ref : (j >= n_kicks ? E1_ref : E0_ref + (j - 1/2) * dE_ref / n_cells)

  dt_ref = (w.outer_in ? (w.outer_in_full ? L_outer : w.outer_in_len)/E_to_v(species, E0_ref) : 0) + 
           (w.outer_out ? (w.outer_out_full ? L_outer : w.outer_out_len)/E_to_v(species, E1_ref) : 0)
  if L_active == 0
    # No inside drifts
  elseif n_cells == 0
    L_inner = L_active / 2
    dt_ref += (w.pre ? (w.pre_full ? L_inner : w.pre_len)/E_to_v(species, E_after(w.ka)) : 0) + 
              (w.post ? (w.post_full ? L_inner : w.post_len)/E_to_v(species, E_after(w.kb+1)) : 0)
  else
    if w.pre
      E_now_ref = E_after(w.ka)
      dt_ref += w.pre_full ? L_active / (n_cells * E_to_v(species, E_now_ref)) : w.pre_len/E_to_v(species, E_now_ref)
    end
    for i_step = w.ka+1:w.kb
      E_now_ref = E0_ref + (i_step - 1/2) * dE_ref / n_cells
      dt_ref += L_active / (n_cells * E_to_v(species, E_now_ref))
    end
    if w.post
      E_now_ref = E_after(w.kb+1)
      dt_ref += w.post_full ? L_active / (n_cells * E_to_v(species, E_now_ref)) : w.post_len/E_to_v(species, E_now_ref)
    end
  end
  return dt_ref, E_after(w.ka), E_after(w.kb+1)
end
