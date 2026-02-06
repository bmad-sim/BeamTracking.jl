# =========== BYPASS UNPACKING FOR RUNGEKUTTA ============= #

"""
Specialized _track! for RungeKutta that bypasses the unpacking system.
Gets field function directly from Beamlines.field_calc and passes the full element.
"""
function _track!(
  coords::Coords,
  bunch::Bunch,
  ele::LineElement,
  tm::RungeKutta,
  scalar_params::Bool,
  ramp_without_rf;
  kwargs...
)
  # Get basic element properties (type-unstable unpacking)
  L = float(ele.L)
  ap = deval(ele.AlignmentParams)
  bp = deval(ele.BendParams)
  dp = deval(ele.ApertureParams)
  patch = deval(ele.PatchParams)
  bm = deval(ele.BMultipoleParams)
  p_over_q_ref = bunch.p_over_q_ref

  if scalar_params
    L = scalarize(L)
    ap = scalarize(ap)
    bp = scalarize(bp)
    bm = scalarize(bm)
    pp = scalarize(pp)
    dp = scalarize(dp)
    mp = scalarize(mp)
    rp = scalarize(rp)
    lp = scalarize(lp)
    p_over_q_ref = scalarize(p_over_q_ref)
  end
  
  # Function barrier
  runge_kutta_universal!(coords, tm, ramp_without_rf, bunch, L, p_over_q_ref, ap, bp, dp, patch, bm; kwargs...)
end

# Step 2: Type-stable computation -----------------------------------------
function runge_kutta_universal!(
  coords,
  tm,
  ramp_without_rf,
  bunch,
  L,
  p_over_q_ref,
  alignmentparams,
  bendparams,
  apertureparams,
  patchparams,
  bmultipoleparams;
  kwargs...
)
  species = bunch.species
  # Setup reference state
  beta_gamma_ref = R_to_beta_gamma(species, p_over_q_ref)
  kc = KernelChain(Val{6}(), RefState(bunch.t_ref, beta_gamma_ref))

  # Evolve time through whole element
  bunch.t_ref += L/beta_gamma_to_v(beta_gamma_ref)

  # Handle reference momentum ramping
  if p_over_q_ref isa TimeDependentParam
    p_over_q_ref_initial = p_over_q_ref
    p_over_q_ref_final = p_over_q_ref(bunch.t_ref)
    if !(p_over_q_ref_initial ≈ p_over_q_ref_final)
      kc = push(kc, KernelCall(BeamTracking.update_P0!, (p_over_q_ref_initial, p_over_q_ref_final, ramp_without_rf)))
      setfield!(bunch, :p_over_q_ref, p_over_q_ref_final)
    end
  end

  # Error conditions
  if L <= 0.0
    error("RungeKutta tracking does not support zero-length elements") 
  end
  if isactive(patchparams)
    error("RungeKutta tracking does not support patch elements")
  end

  # Entrance aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, true)))
        kc = push(kc, @inline(aperture(tm, bunch, apertureparams, true)))
      else
        kc = push(kc, @inline(aperture(tm, bunch, apertureparams, true)))
        kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, true)))
      end
    else
      kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, true)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, bunch, apertureparams, true)))
  end

  # Setup physics parameters
  p_over_q_ref = bunch.p_over_q_ref
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(species, p_over_q_ref)
  charge = chargeof(species)
  p0c = BeamTracking.R_to_pc(species, p_over_q_ref)
  mc2 = massof(species)

  # Determine step size to use
  if tm.ds_step > 0
    ds_step = tm.ds_step
  elseif tm.n_steps > 0
    ds_step = L / tm.n_steps
  else
    ds_step = BeamTracking.DEFAULT_RK4_DS_STEP
  end

  s_span = (0.0, L)

  # Get curvature from BendParams if present
  g_bend = isactive(bendparams) ? bendparams.g_ref : 0.0

  # Extract multipole parameters
  if isactive(bmultipoleparams)
    mm = bmultipoleparams.order
    kn, ks = get_strengths(bmultipoleparams, L, p_over_q_ref)
  else
    # Default to drift
    mm = SVector{0, Int}()
    kn = SVector{0, typeof(L)}()
    ks = SVector{0, typeof(L)}()
  end

  # Build RK4 kernel call
  params = (beta_0, gamsqr_0, tilde_m, charge, p0c, mc2, s_span, ds_step, g_bend, mm, kn, ks)
  kc = push(kc, KernelCall(BeamTracking.RungeKuttaTracking.rk4_kernel!, params))

  # Exit aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = push(kc, @inline(aperture(tm, bunch, apertureparams, false)))
        kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, false)))
      else
        kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, false)))
        kc = push(kc, @inline(aperture(tm, bunch, apertureparams, false)))
      end
    else
      kc = push(kc, @inline(alignment(tm, bunch, alignmentparams, bendparams, L, false)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, bunch, apertureparams, false)))
  end

  # Launch kernels
  @noinline launch!(coords, kc; kwargs...)
  return nothing
end

# =========== ALIGNMENT AND APERTURE ============= #

@inline alignment(tm::RungeKutta, bunch, alignmentparams, bendparams, L, entering) =
  alignment(Exact(), bunch, alignmentparams, bendparams, L, entering)

@inline aperture(tm::RungeKutta, bunch, apertureparams, entering) =
  aperture(Exact(), bunch, apertureparams, entering)
