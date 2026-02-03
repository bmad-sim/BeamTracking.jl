# Step 1: Unpack the element ---------------------------------------------
function _track!(
  coords::Coords,
  bunch::Bunch,
  ele::LineElement, 
  tm,
  scalar_params,
  ramp_without_rf;
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
  p_over_q_ref = lp.beamline.p_over_q_ref

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
  universal!(coords, tm, ramp_without_rf, bunch, L, p_over_q_ref, ap, bp, bm, pp, dp, rp, lp, mp; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  coords,
  tm,
  ramp_without_rf, 
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
  mapparams;
  kwargs...
) 
  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.p_over_q_ref)
  # Current KernelChain length is 6 because we have up to
  # 2 aperture, 2 alignment, 1 body kernel, and 
  # 1 kernel to update the particles' reference energy
  kc = KernelChain(Val{6}(), RefState(bunch.t_ref, beta_gamma_ref))

  # Ramping
  if p_over_q_ref isa TimeDependentParam
    p_over_q_ref_initial = bunch.p_over_q_ref
    p_over_q_ref_final = p_over_q_ref(bunch.t_ref)
    if !(p_over_q_ref_initial ≈ p_over_q_ref_final)
      kc = push(kc, KernelCall(BeamTracking.update_P0!, (p_over_q_ref_initial, p_over_q_ref_final, ramp_without_rf)))
      setfield!(bunch, :p_over_q_ref, p_over_q_ref_final)
    end
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

  if isactive(mapparams)    
    if isactive(alignmentparams)
      error("Tracking through a LineElement containing both MapParams and AlignmentParams is undefined")
    elseif isactive(bendparams)
      error("Tracking through a LineElement containing both MapParams and BendParams not currently defined")
    elseif isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both MapParams and BMultipoleParams not currently defined")
    elseif isactive(rfparams)
      error("Tracking through a LineElement containing both MapParams and RFParams not currently defined")
    elseif isactive(patchparams)
      error("Tracking through a LineElement containing both MapParams and PatchParams not currently defined")
    else
      kc = push(kc, @inline(pure_map(tm, bunch, mapparams, L)))
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
      kc = push(kc, @inline(pure_patch(tm, bunch, patchparams, L)))
    end

  elseif isactive(rfparams)
    !rfparams.is_crabcavity || error("Crab cavities not yet supported for tracking")
    omega = rf_omega(rfparams, beamlineparams.beamline.line[end].s_downstream, bunch.species, bunch.p_over_q_ref)
    t0 = rf_phi0(rfparams) / omega

    if isactive(bendparams)
      error("Tracking through a LineElement containing both RFParams and BendParams not currently defined")
    else
      if !isactive(bmultipoleparams)
        kc = push(kc, @inline(pure_rf(tm, bunch, rfparams, omega, t0, L)))
      else
        kc = push(kc, @inline(bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, omega, t0, L)))
      end
    end

  elseif isactive(bendparams)
    if bendparams.edge1_int != 0 || bendparams.edge2_int != 0; error("edge1_int and edge2_int not yet handled for tracking"); end
    # Bend
    if !isactive(bmultipoleparams) 
      # Bend no field
      kc = push(kc, @inline(bend_no_field(tm, bunch, bendparams, L)))
    else
      n_multipoles = get_n_multipoles(bmultipoleparams)
      if 0 in bmultipoleparams.order # Bend-solenoid
        if n_multipoles == 1
          bm0 = first(bmultipoleparams)
          # Pure bend-solenoid
          kc = push(kc, @inline(bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)))
        else
          # Bend-solenoid with other multipoles of order > 0
          kc = push(kc, @inline(bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)))
        end
      elseif 1 in bmultipoleparams.order # Bend-dipole
        if n_multipoles == 1
          bm1 = first(bmultipoleparams)
          # Pure bend-dipole
          kc = push(kc, @inline(bend_pure_bdipole(tm, bunch, bendparams, bm1, L)))
        else
          # Bend-dipole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)))
        end
      elseif 2 in bmultipoleparams.order # Bend-quadrupole
        if n_multipoles == 1
          bm2 = first(bmultipoleparams)
          # Pure bend-quadrupole
          kc = push(kc, @inline(bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)))
        else
          # Bend-quadrupole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L)))
        end
      else # Bend-multipole
        if n_multipoles == 1
          bmk = first(bmultipoleparams)
          # Pure bend-multipole
          kc = push(kc, @inline(bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)))
        else
          # Bend-multipole with other multipoles of order > 2
          kc = push(kc, @inline(bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)))
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
        kc = push(kc, @inline(pure_bsolenoid(tm, bunch, bm0, L)))
      else
        # Solenoid with other multipoles of order > 0
        kc = push(kc, @inline(bsolenoid(tm, bunch, bmultipoleparams, L)))
      end
    elseif 1 in bmultipoleparams.order # Dipole without bend
      if n_multipoles == 1
        # Pure dipole
        bm1 = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bdipole(tm, bunch, bm1, L)))
      else
        # Dipole with other multipoles of order > 1
        kc = push(kc, @inline(bdipole(tm, bunch, bmultipoleparams, L)))
      end
    elseif 2 in bmultipoleparams.order # Quadrupole
      if n_multipoles == 1
        # Pure quadrupole
        bm2 = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bquadrupole(tm, bunch, bm2, L)))
      else
        # Quadrupole with other multipoles of order > 1
        kc = push(kc, @inline(bquadrupole(tm, bunch, bmultipoleparams, L)))
      end
    else # Higher order multipole
      if n_multipoles == 1
        # Pure multipole
        bmk = first(bmultipoleparams)
        kc = push(kc, @inline(pure_bmultipole(tm, bunch, bmk, L)))
      else
        # Multipole with other multipoles of order > 2
        kc = push(kc, @inline(bmultipole(tm, bunch, bmultipoleparams, L)))
      end
    end

  elseif L != 0
    kc = push(kc, @inline(drift(tm, bunch, L)))
  end

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

  # Evolve time through whole element
  bunch.t_ref += L/beta_gamma_to_v(beta_gamma_ref)
  return kc
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(coords, kc; kwargs...)
  return nothing
end


# === Drift === #
@inline drift(tm, bunch, L) = error("Undefined for tracking method $tm")

# === Coordinate transformations === #
@inline pure_patch(tm, bunch, patchparams, L) = error("Undefined for tracking method $tm")
@inline pure_map(tm, bunch, mapparams, L) = error("Undefined for tracking method $tm")

# === Straight Elements === #
# "Pure" means only ONE SINGLE multipole.
# When "pure" is not present, it means that at least one HIGHER ORDER
# multipole exists.
@inline thin_pure_rf(tm, bunch, rfparams)                          = error("Undefined for tracking method $tm")
@inline thin_pure_bsolenoid(tm, bunch, bm0)                        = error("Undefined for tracking method $tm")
@inline thin_bsolenoid(tm, bunch, bmultipoleparams)                = error("Undefined for tracking method $tm")
@inline thin_pure_bdipole(tm, bunch, bm1)                          = error("Undefined for tracking method $tm")
@inline thin_bdipole(tm, bunch, bmultipoleparams)                  = error("Undefined for tracking method $tm")
@inline thin_pure_bquadrupole(tm, bunch, bm2)                      = error("Undefined for tracking method $tm")
@inline thin_bquadrupole(tm, bunch, bmultipoleparams)              = error("Undefined for tracking method $tm")
@inline thin_pure_bmultipole(tm, bunch, bmk)                       = error("Undefined for tracking method $tm")
@inline thin_bmultipole(tm, bunch, bmultipoleparams)               = error("Undefined for tracking method $tm")
@inline thin_bmultipole_rf(tm, bunch, bmultipoleparams, rfparams)  = error("Undefined for tracking method $tm")

@inline thick_pure_rf(tm, bunch, rfparams, beamlineparams, L)                         = error("Undefined for tracking method $tm")
@inline thick_pure_bsolenoid(tm, bunch, bm0, L)                                       = error("Undefined for tracking method $tm")
@inline thick_bsolenoid(tm, bunch, bmultipoleparams, L)                               = error("Undefined for tracking method $tm")
@inline thick_pure_bdipole(tm, bunch, bm1, L)                                         = error("Undefined for tracking method $tm")
@inline thick_bdipole(tm, bunch, bmultipoleparams, L)                                 = error("Undefined for tracking method $tm")
@inline thick_pure_bquadrupole(tm, bunch, bm2, L)                                     = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, bunch, bmultipoleparams, L)                             = error("Undefined for tracking method $tm")
@inline thick_pure_bmultipole(tm, bunch, bmk, L)                                      = error("Undefined for tracking method $tm")
@inline thick_bmultipole(tm, bunch, bmultipoleparams, L)                              = error("Undefined for tracking method $tm")
@inline thick_bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, beamlineparams, L) = error("Undefined for tracking method $tm")
#=
@inline bend_entrance_fringe(tm, bunch, bendparams, bmultipoleparams, L) = error("Undefined for tracking method $tm")
@inline bend_exit_fringe(tm, bunch, bendparams, bmultipoleparams, L)     = error("Undefined for tracking method $tm")
=#

# === Elements with curving coordinate system "bend" === #
# "Bend" means ONLY a coordinate system curvature through the element.
# It does NOT IMPLY ANY DIPOLE FIELD! Bend specifies if the integration 
# path is curving but does not IMPACT THE PHYSICS inside.
# "Pure" means only ONE SINGLE MULTIPOLE
# When "pure" is not present, it means that at least one higher order 
# multipole exists.

# SciBmad will probably not support thin bends ever but I leave them here for now
@inline thin_bend_no_field(tm, bunch, bendparams)                      = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bsolenoid(tm, bunch, bendparams, bm0)           = error("Undefined for tracking method $tm")
@inline thin_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams)   = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bdipole(tm, bunch, bendparams, bm1)             = error("Undefined for tracking method $tm")
@inline thin_bend_bdipole(tm, bunch, bendparams, bmultipoleparams)     = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bquadrupole(tm, bunch, bendparams, bm2)         = error("Undefined for tracking method $tm")
@inline thin_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams) = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bmultipole(tm, bunch, bendparams, bmk)          = error("Undefined for tracking method $tm")
@inline thin_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams)  = error("Undefined for tracking method $tm")

@inline thick_bend_no_field(tm, bunch, bendparams, L)                      = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)           = error("Undefined for tracking method $tm")
@inline thick_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)   = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bdipole(tm, bunch, bendparams, bm1, L)             = error("Undefined for tracking method $tm")
@inline thick_bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)     = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)         = error("Undefined for tracking method $tm")
@inline thick_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L) = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)          = error("Undefined for tracking method $tm")
@inline thick_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)  = error("Undefined for tracking method $tm")


# === Elements thin vs thick check === #
@inline pure_rf(tm, bunch, rfparams, omega, t0, L)                          = L == 0 ? thin_pure_rf(tm, bunch, rfparams, omega, t0)                         : thick_pure_rf(tm, bunch, rfparams, omega, t0, L)
@inline pure_bsolenoid(tm, bunch, bm0, L)                                   = L == 0 ? thin_pure_bsolenoid(tm, bunch, bm0)                                  : thick_pure_bsolenoid(tm, bunch, bm0, L)      
@inline bsolenoid(tm, bunch, bmultipoleparams, L)                           = L == 0 ? thin_bsolenoid(tm, bunch, bmultipoleparams)                          : thick_bsolenoid(tm, bunch, bmultipoleparams, L)       
@inline pure_bdipole(tm, bunch, bm1, L)                                     = L == 0 ? thin_pure_bdipole(tm, bunch, bm1)                                    : thick_pure_bdipole(tm, bunch, bm1, L)          
@inline bdipole(tm, bunch, bmultipoleparams, L)                             = L == 0 ? thin_bdipole(tm, bunch, bmultipoleparams)                            : thick_bdipole(tm, bunch, bmultipoleparams, L)             
@inline pure_bquadrupole(tm, bunch, bm2, L)                                 = L == 0 ? thin_pure_bquadrupole(tm, bunch, bm2)                                : thick_pure_bquadrupole(tm, bunch, bm2, L)        
@inline bquadrupole(tm, bunch, bmultipoleparams, L)                         = L == 0 ? thin_bquadrupole(tm, bunch, bmultipoleparams)                        : thick_bquadrupole(tm, bunch, bmultipoleparams, L)           
@inline pure_bmultipole(tm, bunch, bmk, L)                                  = L == 0 ? thin_pure_bmultipole(tm, bunch, bmk)                                 : thick_pure_bmultipole(tm, bunch, bmk, L)                   
@inline bmultipole(tm, bunch, bmultipoleparams, L)                          = L == 0 ? thin_bmultipole(tm, bunch, bmultipoleparams)                         : thick_bmultipole(tm, bunch, bmultipoleparams, L)
@inline bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, omega, t0, L)  = L == 0 ? thin_bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, omega, t0) : thick_bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, omega, t0, L)                           
@inline bend_no_field(tm, bunch, bendparams, L)                             = L == 0 ? thin_bend_no_field(tm, bunch, bendparams)                            : thick_bend_no_field(tm, bunch, bendparams, L)
@inline bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)                  = L == 0 ? thin_bend_pure_bsolenoid(tm, bunch, bendparams, bm0)                 : thick_bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)      
@inline bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)          = L == 0 ? thin_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams)         : thick_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)         
@inline bend_pure_bdipole(tm, bunch, bendparams, bm1, L)                    = L == 0 ? thin_bend_pure_bdipole(tm, bunch, bendparams, bm1)                   : thick_bend_pure_bdipole(tm, bunch, bendparams, bm1, L)          
@inline bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)            = L == 0 ? thin_bend_bdipole(tm, bunch, bendparams, bmultipoleparams)           : thick_bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)             
@inline bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)                = L == 0 ? thin_bend_pure_bquadrupole(tm, bunch, bendparams, bm2)               : thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)        
@inline bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L)        = L == 0 ? thin_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams)       : thick_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L)           
@inline bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)                 = L == 0 ? thin_bend_pure_bmultipole(tm, bunch, bendparams, bmk)                : thick_bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)                   
@inline bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)         = L == 0 ? thin_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams)        : thick_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)                      