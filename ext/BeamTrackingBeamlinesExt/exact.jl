 # Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  v,
  work,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::Exact;
  kwargs...
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  pp = ele.PatchParams
  L = ele.L

  # Function barrier
  exact_universal!(i, v, work, bunch, L, bm, bp, ma, pp; kwargs...)
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
  alignmentparams,
  patchparams;
  kwargs...
) 
  
  if isactive(patchparams) # Patch
    if isactive(bendparams) || isactive(bmultipoleparams)
      error("Patch cannot yet be used with bend or multipole elements")
    else
      tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
      winv = ExactTracking.w_inv_matrix(patchparams.dx_rot, patchparams.dy_rot, patchparams.dz_rot)
      runkernel!(ExactTracking.patch!, i, v, work, tilde_m, patchparams.dt, patchparams.dx, patchparams.dy, patchparams.dz, winv, L)
    end
  elseif isactive(bendparams) # Bend
    if !isactive(bmultipoleparams) # Exact bend
      error("Exact bend not yet implemented")
    else # Combined function bend
      error("Exact bend cannot be used with multipole elements")
    end
  elseif isactive(bmultipoleparams) # Straight multipole
    if haskey(bmultipoleparams.bdict, 0) # Solenoid
      if L == 0
        error("Exact thin-lens solenoid not yet implemented (L = 0)")
      else
        Ks = get_thick_strength(bmultipoleparams.bdict[0], L, bunch.Brho_ref)
        tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
        runkernel!(ExactTracking.exact_solenoid!, i, v, work, Ks, beta_0, gamsqr_0, tilde_m, L)
      end
    elseif haskey(bmultipoleparams.bdict, 1) # Kick
      error("Exact kick not yet implemented")
    else # Multipole
      error("Exact solution does not exist for multipole elements")
    end
  elseif L != 0 # Drift
    tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.Brho_ref)
    runkernel!(ExactTracking.exact_drift!, i, v, work, beta_0, gamsqr_0, tilde_m, L)
  end

  return v
end
  
