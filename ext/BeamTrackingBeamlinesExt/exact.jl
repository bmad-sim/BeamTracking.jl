
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
  
  # Step 1: Unpack the element ---------------------------------------------
  function _track!(
    i,
    v,
    work,
    bunch::Bunch,
    ele::Union{LineElement,BitsLineElement}, 
    ::Exact
  )
    # Unpack the line element
    bm = ele.BMultipoleParams
    bp = ele.BendParams
    pp = ele.PatchParams
    L = ele.L
  
    # Function barrier
    exact_universal!(i, v, work, bunch, L, bm, bp, pp)
  end
  
  
  # Step 2: Push particles through -----------------------------------------
  function exact_universal!(
    i, 
    v, 
    work,
    bunch,
    L,  
    bmultipoleparams,
    bendparams,
    patchparams
  ) 
    
    if isactive(patchparams) # Patch
      if isactive(bendparams) || isactive(bmultipoleparams)
        error("Patch cannot yet be used with bend or multipole elements")
      else
        p0c = bunch.Brho_ref*C_LIGHT*bunch.charge/Q
        winv = w_inv_matrix(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
        runkernel!(ExactTracking.patch!, i, v, work, L, p0c, massof(bunch.species), patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv)
      end
    elseif isactive(bendparams) # Bend
      if isactive(bmultipoleparams) # Combined function bend
        error("Exact bend cannot be used with multipole elements")
      else
      end
    elseif isactive(bmultipoleparams) # Straight multipole
      if haskey(bmultipoleparams.bdict, 0) # Solenoid
        if L < 1e-6
          error("Exact thin-lens solenoid not yet implemented (L = 0)")
        else
          Ks = get_thick_strength(bmultipole.bdict[0], L, bunch.Brho_ref)*chargeof(bunch.species)
          p0c = bunch.Brho_ref*C_LIGHT*chargeof(bunch.species)
          tilde_m = massof(bunch.species)/p0c
          gamsqr_0 = 1 + 1/tilde_m^2
          beta_0 = sqrt(1 - 1/gamsqr_0)
          runkernel!(ExactTracking.exact_solenoid!, i, v, work, L, Ks, tilde_m, gamsqr_0, beta_0)
        end
      elseif haskey(bmultipoleparams.bdict, 1) # Kick
        error("Exact kick not yet implemented")
      else # Multipole
        error("Exact solution does not exist for multipole elements")
      end
    elseif L != 0 # Drift
      p0c = bunch.Brho_ref*C_LIGHT*chargeof(bunch.species)
      tilde_m = massof(bunch.species)/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = sqrt(1 - 1/gamsqr_0)
      runkernel!(ExactTracking.exact_drift!, i, v, work, L, tilde_m, gamsqr_0, beta_0)
    end
  
    return v
  end
  