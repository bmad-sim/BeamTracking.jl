
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]


# =========== STRAIGHT ELEMENTS ============= #
# === Thin elements === #
# Only thin element currently supported is thin quadrupole
@inline function thin_pure_bdipole(tm::Linear, bunch, bm1)
  # This a thin (0 length) corrector coil
  # In Fortran Bmad is g == 0 but dgL != 0. 
  # In SciBmad this is g == 0 but Kn0L != 0.
  error("Undefined for tracking method $tm")
end

@inline function thin_bdipole(tm::Linear, bunch, bmultipoleparams)
  if 2 in bmultipoleparams.order
    # This is a thin corrector coil with a quadrupole term
    # In Fortran Bmad is g == 0 but dgL != 0 and Kn1L != 0
    # In SciBmad this is g == 0 but Kn0L != 0 and Kn1L != 0
    error("Undefined for tracking method $tm")
  else # ignore higher order multipoles
    return thin_pure_bdipole(tm, bunch, bmultipoleparams[1])
  end
end

@inline function thin_pure_bquadrupole(tm::Linear, bunch, bm2)
  Kn1L, Ks1L = get_integrated_strengths(bm2, 0, bunch.R_ref)
  # Temporary
  Ks1L ≈ 0 || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
  mx, my = LinearTracking.linear_thin_quad_matrices(Kn1L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, 0, nothing, nothing))
end

# Ignore higher order multipoles with quadrupole: 
@inline thin_bquadrupole(tm::Linear, bunch, bmultipoleparams) = thin_pure_bquadrupole(tm, bunch, bmultipoleparams[2])

# Ignore higher order multipoles: do nothing for thin case (0 length)
@inline thin_pure_bmultipole(tm::Linear, bunch, bmk) = KernelCall()
@inline thin_bmultipole(tm::Linear, bunch, bmultipoleparams) = KernelCall()


# === Thick elements === #
@inline function drift(tm::Linear, bunch, L)
  gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
  return KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2))
end

@inline function thick_pure_bsolenoid(tm::Linear, bunch, bm0, L)
  gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
  Ksol, _ = get_strengths(bm0, L, bunch.R_ref)
  mxy = LinearTracking.linear_solenoid_matrix(Ksol, L)
  return KernelCall(LinearTracking.linear_coast!, (mxy, L/gamma_0^2, nothing, nothing))
end

@inline function thick_pure_bdipole(tm::Linear, bunch, bm1, L)
  # This a thick corrector coil
  # In Fortran Bmad is g == 0 but dg != 0. 
  # In SciBmad this is g == 0 but Kn0 != 0.
  # Singularity for g==0
  error("Undefined for tracking method $tm")
end

@inline function thick_bdipole(tm::Linear, bunch, bmultipoleparams, L)
  if 2 in bmultipoleparams.order
    # This is a thick corrector coil with a quadrupole term
    # In Fortran Bmad is g == 0 but dg != 0 and Kn1 != 0
    # In SciBmad this is g == 0 but Kn0 != 0 and Kn1 != 0
    gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
    Kn0, Ks0 = get_strengths(bmultipoleparams[1], L, bunch.R_ref)
    Kn1, Ks1 = get_strengths(bmultipoleparams[2], L, bunch.R_ref) 
    (Ks0 ≈ 0 && Ks1 ≈ 0) || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
    mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(0, 0, 0, Kn0, Kn1, gamma_0, L)
    return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
  else # ignore higher order multipoles
    return thick_pure_bdipole(tm, bunch, bmultipoleparams[1], L)
  end
end

@inline function thick_pure_bquadrupole(tm::Linear, bunch, bm2, L)
  gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
  Kn1, Ks1 = get_strengths(bm2, L, bunch.R_ref)
  Ks1 ≈ 0 || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
  mx, my = LinearTracking.linear_quad_matrices(Kn1, L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, L/gamma_0^2, nothing, nothing))
end

# Ignore higher order multipoles with quadrupole:
@inline thick_bquadrupole(tm::Linear, bunch, bmultipoleparams, L) = thick_pure_bquadrupole(tm, bunch, bmultipoleparams[2], L)

# Ignore higher order multipoles: drift
@inline thick_pure_bmultipole(tm::Linear, bunch, bmk, L) = drift(tm, bunch, L)
@inline thick_bmultipole(tm::Linear, bunch, bmultipoleparams, L) = drift(tm, bunch, L)

# =========== BENDING ELEMENTS ============= #
# Again, "bend" means coordinate system bends, nothing wrt the field.
# SciBmad likely will not support "thin" bends

@inline function thick_bend_no_field(tm::Linear, bunch, bendparams, L)
  # This has NO FIELD!
  # In Fortran Bmad it is like setting dg == -g
  # Singularity for simutaneous`Kn0==0' and 'Kn1 = 0'
  error("Undefined for tracking method $tm")
end

@inline function thick_bend_pure_bdipole(tm::Linear, bunch, bendparams, bm1, L)     
  gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
  Kn0, Ks0 = get_strengths(bm1, L, bunch.R_ref)
  Ks0 ≈ 0 || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
  mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(bendparams.g_ref, bendparams.e1, bendparams.e2, Kn0, 0, gamma_0, L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
end

@inline function thick_bend_bdipole(tm::Linear, bunch, bendparams, bmultipoleparams, L)   
  if 2 in bmultipoleparams.order
    # This is a thick combined function magnet
    # In Fortran Bmad is g != 0, dg != 0, Kn1 != 0
    # In SciBmad this is g != 0, Kn0 != 0, Kn1 != 0
    gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
    Kn0, Ks0 = get_strengths(bmultipoleparams[1], L, bunch.R_ref)
    Kn1, Ks1 = get_strengths(bmultipoleparams[2], L, bunch.R_ref)
    (Ks0 ≈ 0 && Ks1 ≈ 0) || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
    mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(bendparams.g_ref, bendparams.e1, bendparams.e2, Kn0, Kn1, gamma_0, L)
    return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
  else # ignore higher order multipoles
    return thick_bend_pure_bdipole(tm, bunch, bendparams, bmultipoleparams[1], L)
  end
end

@inline function thick_bend_pure_bquadrupole(tm::Linear, bunch, bendparams, bm2, L) 
  # This is a quadrupole with a g, I think your code should be able to handle this
  # In Fortran Bmad it would be like dg == -g, Kn1 != 0.
  gamma_0 = R_to_gamma(bunch.species, bunch.R_ref)
  Kn1, Ks1 = get_strengths(bm2, L, bunch.R_ref)
  Ks1 ≈ 0 || error("Skew/tilt multipoles not implemented yet for tracking method $tm")
  mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(bendparams.g_ref, bendparams.e1, bendparams.e2, 0, Kn1, gamma_0, L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
end

# Ignore higher order multipoles:
@inline thick_bend_bquadrupole(tm::Linear, bunch, bendparams, bmultipoleparams, L) = thick_bend_pure_bquadrupole(tm, bunch, bendparams, bmultipoleparams[2], L)
# Ignore higher order multipoles: treat like bend no field
@inline thick_bend_pure_bmultipole(tm::Linear, bunch, bendparams, bmk, L) = thick_bend_no_field(tm, bunch, bendparams, L)
@inline thick_bend_bmultipole(tm::Linear, bunch, bendparams, bmultipoleparams, L) = thick_bend_no_field(tm, bunch, bendparams, L)
