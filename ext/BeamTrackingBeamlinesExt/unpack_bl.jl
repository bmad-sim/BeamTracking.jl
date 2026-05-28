# Step 1: Unpack the element ---------------------------------------------
function _track!(
  coords::Coords,
  bunch::Bunch,
  ele::LineElement, 
  p_over_q_ref,
  tm,
  scalar_params,
  ramp_particle_energy_without_rf,
  ramp_update_each_particle;
  kwargs...
)
  # Unpack the line element (type unstable)
  L = float(ele.L) # Automatically calls deval (element-level get)
  # float call is required because L is allowed to be any type
  # in order to keep binaries smaller for tracking routines, 
  # we don't want to compile separate routines for Int64
  ap = deval(ele.AlignmentParams)
  bp = deval(ele.BendParams)
  bm = deval(ele.BMultipoleParams)
  pp = deval(ele.PatchParams)
  dp = deval(ele.ApertureParams)
  mp = deval(ele.MapParams)
  rp = deval(ele.RFParams)
  lp = deval(ele.BeamlineParams)
  fpp = deval(ele.FourPotentialParams)

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
    fpp = scalarize(fpp)
    p_over_q_ref = scalarize(p_over_q_ref)
  end

  # Function barrier
  universal!(coords, tm, ele, ramp_particle_energy_without_rf, ramp_update_each_particle, bunch, L, p_over_q_ref, ap, bp, bm, pp, dp, rp, lp, mp, fpp; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  coords,
  tm,
  ele,
  ramp_particle_energy_without_rf, 
  ramp_update_each_particle,
  bunch,
  L, 
  p_over_q_ref,
  alignmentparams,
  bendparams,
  bmultipoleparams,
  patchparams,
  apertureparams,
  rfparams,
  beamlineparams,
  mapparams,
  fourpotentialparams;
  kwargs...
) 
  if p_over_q_ref isa TimeDependentParam
    beta_gamma_ref0 = R_to_beta_gamma(bunch.species, p_over_q_ref(bunch.t_ref))
  else
    beta_gamma_ref0 = R_to_beta_gamma(bunch.species, p_over_q_ref)
  end

  # Current KernelChain length is 9 because we have up to
  # 2 aperture, 2 alignment, 1 body kernel, 1 IBS kernel,
  # 2 kernels to update the particles' reference energy,
  # and 2 for coordinate conversion with implicit
  kc = KernelChain(Val{10}(), RefState(bunch.t_ref, beta_gamma_ref0))

  p_over_q_ref_initial = bunch.p_over_q_ref
  ramp_per_particle = p_over_q_ref isa TimeDependentParam && ramp_update_each_particle

  if ramp_per_particle
    kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch.p_over_q_ref, p_over_q_ref-bunch.p_over_q_ref, Val{!ramp_particle_energy_without_rf}())))
  else
    # Make sure to evaluate p_over_q_ref if not ramp_update_each_particle
    p_over_q_ref = p_over_q_ref isa TimeDependentParam ? p_over_q_ref(bunch.t_ref) : p_over_q_ref
    if !(p_over_q_ref ≈ p_over_q_ref_initial)
        kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (p_over_q_ref_initial, p_over_q_ref - p_over_q_ref_initial, Val{!ramp_particle_energy_without_rf}())))
        bunch.p_over_q_ref = p_over_q_ref
    end
  end

  # Entrance aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
      else
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
      end
    else
      kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
  end

  if ((hasfield(typeof(tm), :ibs_damping_on) && hasfield(typeof(tm), :ibs_fluctuations_on)) 
    && (tm.ibs_damping_on || tm.ibs_fluctuations_on) && L > 0)
    bp = ifelse(isactive(bendparams), bendparams, nothing)
    kc = push(kc, @inline(ibs_kick(tm, p_over_q_ref, bunch, bp, L)))
  end

  if isactive(mapparams)    
    if isactive(bendparams)
      error("Tracking through a LineElement containing both MapParams and BendParams not currently defined")
    elseif isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both MapParams and BMultipoleParams not currently defined")
    elseif isactive(rfparams)
      error("Tracking through a LineElement containing both MapParams and RFParams not currently defined")
    elseif isactive(patchparams)
      error("Tracking through a LineElement containing both MapParams and PatchParams not currently defined")
    elseif isactive(fourpotentialparams)
      error("Tracking through a LineElement containing both MapParams and FourPotentialParams not currently defined")
    else
      kc = push(kc, @inline(pure_map(tm, p_over_q_ref, bunch, mapparams, L)))
    end

  elseif isactive(fourpotentialparams)    
    if isactive(alignmentparams)
      error("Tracking through a LineElement containing both FourPotentialParams and AlignmentParams is undefined")
    elseif isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both FourPotentialParams and BMultipoleParams not currently defined")
    elseif isactive(rfparams)
      error("Tracking through a LineElement containing both FourPotentialParams and RFParams not currently defined")
    elseif isactive(patchparams)
      error("Tracking through a LineElement containing both MapParams and PatchParams not currently defined")
    else
      kc = push(kc, @inline(implicit_in(tm, p_over_q_ref, bunch)))
      kc = push(kc, @inline(implicit_body(tm, p_over_q_ref, bunch, fourpotentialparams, bendparams, L)))
      kc = push(kc, @inline(implicit_out(tm, p_over_q_ref, bunch)))
    end

  elseif isactive(patchparams)    
    if isactive(alignmentparams)
      error("Tracking through a LineElement containing both PatchParams and AlignmentParams is undefined")
    elseif isactive(bendparams)
      error("Tracking through a LineElement containing both PatchParams and BendParams not currently defined")
    elseif isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both PatchParams and BMultipoleParams not currently defined")
    elseif isactive(rfparams)
      error("Tracking through a LineElement containing both PatchParams and RFParams not currently defined")
    else
      # Pure patch
      kc = push(kc, @inline(pure_patch(tm, p_over_q_ref, bunch, patchparams, L)))
    end

  elseif isactive(rfparams)
    if isactive(bendparams)
      error("Tracking through a LineElement containing both RFParams and BendParams not currently defined")
    end
    !rfparams.is_crabcavity || error("Crab cavities not yet supported for tracking")

    kc = push(kc, @inline(rfcavity(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L)))
    
  elseif isactive(bendparams)
    if bendparams.edge1_int != 0 || bendparams.edge2_int != 0; error("edge1_int and edge2_int not yet handled for tracking"); end
    # Bend
    if !isactive(bmultipoleparams) 
      # Bend no field
      kc = push(kc, @inline(bend_no_field(tm, p_over_q_ref, bunch, bendparams, L)))
    else
      n_multipoles = get_n_multipoles(bmultipoleparams)
      if 0 in bmultipoleparams.order # Bend-solenoid
        if n_multipoles == 1
          bm0 = first(bmultipoleparams)
          # Pure bend-solenoid
          kc = push(kc, @inline(bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0, L)))
        else
          # Bend-solenoid with other multipoles of order > 0
          kc = push(kc, @inline(bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)))
        end
      elseif 1 in bmultipoleparams.order # Bend-dipole
        if n_multipoles == 1
          bm1 = first(bmultipoleparams)
          # Pure bend-dipole
          kc = push(kc, @inline(bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1, L)))
        else
          # Bend-dipole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)))
        end
      elseif 2 in bmultipoleparams.order # Bend-quadrupole
        if n_multipoles == 1
          bm2 = first(bmultipoleparams)
          # Pure bend-quadrupole
          kc = push(kc, @inline(bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2, L)))
        else
          # Bend-quadrupole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)))
        end
      else # Bend-multipole
        if n_multipoles == 1
          bmk = first(bmultipoleparams)
          # Pure bend-multipole
          kc = push(kc, @inline(bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk, L)))
        else
          # Bend-multipole with other multipoles of order > 2
          kc = push(kc, @inline(bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)))
        end
      end
    end

  elseif isactive(bmultipoleparams)
    # BMultipole
    n_multipoles = get_n_multipoles(bmultipoleparams)
    if 0 in bmultipoleparams.order # Solenoid
      if n_multipoles == 1
        # Pure solenoid
        bm0 = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bsolenoid(tm, p_over_q_ref, bunch, bm0, L)))
      else
        # Solenoid with other multipoles of order > 0
        kc = push(kc, @inline(bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams, L)))
      end
    elseif 1 in bmultipoleparams.order # Dipole without bend
      if n_multipoles == 1
        # Pure dipole
        bm1 = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bdipole(tm, p_over_q_ref, bunch, bm1, L)))
      else
        # Dipole with other multipoles of order > 1
        kc = push(kc, @inline(bdipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)))
      end
    elseif 2 in bmultipoleparams.order # Quadrupole
      if n_multipoles == 1
        # Pure quadrupole
        bm2 = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bquadrupole(tm, p_over_q_ref, bunch, bm2, L)))
      else
        # Quadrupole with other multipoles of order > 1
        kc = push(kc, @inline(bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams, L)))
      end
    else # Higher order multipole
      if n_multipoles == 1
        # Pure multipole
        bmk = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bmultipole(tm, p_over_q_ref, bunch, bmk, L)))
      else
        # Multipole with other multipoles of order > 2
        kc = push(kc, @inline(bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)))
      end
    end

  elseif L != 0
    kc = push(kc, @inline(drift(tm, p_over_q_ref, bunch, L)))
  end

  # Exit aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
      else
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
      end
    else
      kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
  end

  # noinline necessary here for small binaries and faster execution
  @noinline launch!(coords, kc; kwargs...)

  # Reference time evolution thru element assumes constant energy
  # using the energy at the start of the element:
  bunch.t_ref += L / beta_gamma_to_v(beta_gamma_ref0)

  if first(kc.chain).kernel == BeamTracking.blank_kernel! # Still execute callbacks if nothing happened.
    BeamTracking.execute_callbacks(bunch.coords, 0, (0, 0))
  end

  # If ramping, now need to uniformly ramp all particles to same reference energy
  if ramp_per_particle
    # If TimeDependentParam, at end ramp all 
    # uniformly to p_over_q_ref(t_ref at end)
    bunch.p_over_q_ref = p_over_q_ref(bunch.t_ref)
    kc = push(kc, make_kernel_call(BeamTracking.reference_momentum_shift!, (bunch.p_over_q_ref, bunch.p_over_q_ref-p_over_q_ref, Val{!ramp_particle_energy_without_rf}())))
  end

  return nothing
end

#---------------------------------------------------------------------------------------------------
# universal! for SaganCavity tracking.

function universal!(coords, tm::SaganCavity, ele, ramp_particle_energy_without_rf, ramp_update_each_particle, bunch, L,
  p_over_q_ref, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams,
  rfparams, beamlineparams, mapparams, fourpotentialparams; kwargs...) 

  !isactive(mapparams) || error("SaganCavity Tracking through element $ele_name with MapParams is undefined")
  !isactive(patchparams) || error("SaganCavity Tracking through element $ele_name with PatchParams is undefined")
  !isactive(patchparams) || error("SaganCavity Tracking through element $ele_name with BendParams is undefined")
  !isactive(fourpotentialparams) || error("SaganCavity Tracking through element $ele_name with FourPotentialParams is undefined")
  isactive(rfparams) || error("SaganCavity Tracking through element $ele_name without RFParams is undefined")

  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)
  kc = KernelChain(Val{6}(), RefState(bunch.t_ref, beta_gamma_ref))

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
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
      else
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
      end
    else
      kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, true)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, true)))
  end

  # Cavity tracking
  kc = push(kc, @inline(sagan_cavity(tm, p_over_q_ref, bunch, ele.name, bmultipoleparams, rfparams, beamlineparams, L)))

  # Exit aperture and alignment
  if isactive(alignmentparams)
    if isactive(apertureparams)
      if apertureparams.aperture_shifts_with_body
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
      else
        kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
        kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
      end
    else
      kc = push(kc, @inline(alignment(tm, p_over_q_ref, bunch, alignmentparams, bendparams, L, false)))
    end
  elseif isactive(apertureparams)
    kc = push(kc, @inline(aperture(tm, p_over_q_ref, bunch, apertureparams, false)))
  end

  # noinline necessary here for small binaries and faster execution
  @noinline launch!(coords, kc; kwargs...)

  # reference time change
  if L != 0
    species = bunch.species
    p1_over_q_ref = beamlineparams.beamline.p_over_q_ref
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

    bunch.t_ref += dt_ref
  end

  return nothing

end

#---------------------------------------------------------------------------------------------------

# === Drift === #
@inline drift(tm, p_over_q_ref, bunch, L) = error("Undefined for tracking method $tm")

# == Implicit === #
@inline implicit_in(tm, p_over_q_ref, bunch) = error("Undefined for tracking method $tm")
@inline implicit_body(tm, p_over_q_ref, bunch, fourpotentialparams, bendparams, L) = error("Undefined for tracking method $tm")
@inline implicit_out(tm, p_over_q_ref, bunch) = error("Undefined for tracking method $tm")

# === Straight Elements === #
# "Pure" means only ONE SINGLE multipole.
# When "pure" is not present, it means that at least one HIGHER ORDER
# multipole exists.
@inline thin_pure_rf(tm, p_over_q_ref, bunch, rfparams)                          = error("Undefined for tracking method $tm")
@inline thin_pure_bsolenoid(tm, p_over_q_ref, bunch, bm0)                        = error("Undefined for tracking method $tm")
@inline thin_bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams)                = error("Undefined for tracking method $tm")
@inline thin_pure_bdipole(tm, p_over_q_ref, bunch, bm1)                          = error("Undefined for tracking method $tm")
@inline thin_bdipole(tm, p_over_q_ref, bunch, bmultipoleparams)                  = error("Undefined for tracking method $tm")
@inline thin_pure_bquadrupole(tm, p_over_q_ref, bunch, bm2)                      = error("Undefined for tracking method $tm")
@inline thin_bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams)              = error("Undefined for tracking method $tm")
@inline thin_pure_bmultipole(tm, p_over_q_ref, bunch, bmk)                       = error("Undefined for tracking method $tm")
@inline thin_bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams)               = error("Undefined for tracking method $tm")
@inline thin_bmultipole_rf(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams)  = error("Undefined for tracking method $tm")

@inline thick_pure_rf(tm, p_over_q_ref, bunch, rfparams, beamlineparams, L)                         = error("Undefined for tracking method $tm")
@inline thick_pure_bsolenoid(tm, p_over_q_ref, bunch, bm0, L)                                       = error("Undefined for tracking method $tm")
@inline thick_bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams, L)                               = error("Undefined for tracking method $tm")
@inline thick_pure_bdipole(tm, p_over_q_ref, bunch, bm1, L)                                         = error("Undefined for tracking method $tm")
@inline thick_bdipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                                 = error("Undefined for tracking method $tm")
@inline thick_pure_bquadrupole(tm, p_over_q_ref, bunch, bm2, L)                                     = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                             = error("Undefined for tracking method $tm")
@inline thick_pure_bmultipole(tm, p_over_q_ref, bunch, bmk, L)                                      = error("Undefined for tracking method $tm")
@inline thick_bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                              = error("Undefined for tracking method $tm")
@inline thick_bmultipole_rf(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L) = error("Undefined for tracking method $tm")

# === Elements with curving coordinate system "bend" === #
# "Bend" means ONLY a coordinate system curvature through the element.
# It does NOT IMPLY ANY DIPOLE FIELD! Bend specifies if the integration 
# path is curving but does not IMPACT THE PHYSICS inside.
# "Pure" means only ONE SINGLE MULTIPOLE
# When "pure" is not present, it means that at least one higher order 
# multipole exists.

# SciBmad will probably not support thin bends ever but I leave them here for now
@inline thin_bend_no_field(tm, p_over_q_ref, bunch, bendparams)                      = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0)           = error("Undefined for tracking method $tm")
@inline thin_bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)   = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1)             = error("Undefined for tracking method $tm")
@inline thin_bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)     = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2)         = error("Undefined for tracking method $tm")
@inline thin_bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams) = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk)          = error("Undefined for tracking method $tm")
@inline thin_bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)  = error("Undefined for tracking method $tm")

@inline thick_bend_no_field(tm, p_over_q_ref, bunch, bendparams, L)                      = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0, L)           = error("Undefined for tracking method $tm")
@inline thick_bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)   = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1, L)             = error("Undefined for tracking method $tm")
@inline thick_bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)     = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2, L)         = error("Undefined for tracking method $tm")
@inline thick_bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L) = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk, L)          = error("Undefined for tracking method $tm")
@inline thick_bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)  = error("Undefined for tracking method $tm")


# === Elements thin vs thick check === #
@inline pure_rf(tm, p_over_q_ref, bunch, rfparams, beamlineparams, L)                          = L == 0 ? thin_pure_rf(tm, p_over_q_ref, bunch, rfparams, beamlineparams)                         : thick_pure_rf(tm, p_over_q_ref, bunch, rfparams, beamlineparams, L)
@inline pure_bsolenoid(tm, p_over_q_ref, bunch, bm0, L)                                   = L == 0 ? thin_pure_bsolenoid(tm, p_over_q_ref, bunch, bm0)                                  : thick_pure_bsolenoid(tm, p_over_q_ref, bunch, bm0, L)      
@inline bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams, L)                           = L == 0 ? thin_bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams)                          : thick_bsolenoid(tm, p_over_q_ref, bunch, bmultipoleparams, L)       
@inline pure_bdipole(tm, p_over_q_ref, bunch, bm1, L)                                     = L == 0 ? thin_pure_bdipole(tm, p_over_q_ref, bunch, bm1)                                    : thick_pure_bdipole(tm, p_over_q_ref, bunch, bm1, L)          
@inline bdipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                             = L == 0 ? thin_bdipole(tm, p_over_q_ref, bunch, bmultipoleparams)                            : thick_bdipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)             
@inline pure_bquadrupole(tm, p_over_q_ref, bunch, bm2, L)                                 = L == 0 ? thin_pure_bquadrupole(tm, p_over_q_ref, bunch, bm2)                                : thick_pure_bquadrupole(tm, p_over_q_ref, bunch, bm2, L)        
@inline bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                         = L == 0 ? thin_bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams)                        : thick_bquadrupole(tm, p_over_q_ref, bunch, bmultipoleparams, L)           
@inline pure_bmultipole(tm, p_over_q_ref, bunch, bmk, L)                                  = L == 0 ? thin_pure_bmultipole(tm, p_over_q_ref, bunch, bmk)                                 : thick_pure_bmultipole(tm, p_over_q_ref, bunch, bmk, L)                   
@inline bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)                          = L == 0 ? thin_bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams)                         : thick_bmultipole(tm, p_over_q_ref, bunch, bmultipoleparams, L)
@inline bmultipole_rf(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L)  = L == 0 ? thin_bmultipole_rf(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams) : thick_bmultipole_rf(tm, p_over_q_ref, bunch, bmultipoleparams, rfparams, beamlineparams, L)        
@inline bend_no_field(tm, p_over_q_ref, bunch, bendparams, L)                             = L == 0 ? thin_bend_no_field(tm, p_over_q_ref, bunch, bendparams)                            : thick_bend_no_field(tm, p_over_q_ref, bunch, bendparams, L)
@inline bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0, L)                  = L == 0 ? thin_bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0)                 : thick_bend_pure_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bm0, L)      
@inline bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)          = L == 0 ? thin_bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)         : thick_bend_bsolenoid(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)         
@inline bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1, L)                    = L == 0 ? thin_bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1)                   : thick_bend_pure_bdipole(tm, p_over_q_ref, bunch, bendparams, bm1, L)          
@inline bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)            = L == 0 ? thin_bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)           : thick_bend_bdipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)             
@inline bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2, L)                = L == 0 ? thin_bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2)               : thick_bend_pure_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bm2, L)        
@inline bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)        = L == 0 ? thin_bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)       : thick_bend_bquadrupole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)           
@inline bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk, L)                 = L == 0 ? thin_bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk)                : thick_bend_pure_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmk, L)                   
@inline bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)         = L == 0 ? thin_bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams)        : thick_bend_bmultipole(tm, p_over_q_ref, bunch, bendparams, bmultipoleparams, L)                      