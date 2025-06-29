
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]


# =========== STRAIGHT ELEMENTS ============= #
# === Thin elements === #
# Only thin element currently supported is thin quadrupole
@inline function thin_pure_bdipole(tm::Linear, bunch, bm1)
  # Sophia: this a thin (0 length) corrector coil
  # In Fortran Bmad is g == 0 but dgL != 0. 
  # In SciBmad this is g == 0 but K0L != 0.
  error("Undefined for tracking method $tm")
end

@inline function thin_bdipole(tm::Linear, bunch, bdict)
  if haskey(bdict, 2)
    # Sophia: this is a thin corrector coil with a quadrupole term
    # In Fortran Bmad is g == 0 but dgL != 0 and K1L != 0
    # In SciBmad this is g == 0 but K0L != 0 and K1L != 0
    error("Undefined for tracking method $tm")
  else # ignore higher order multipoles
    return thin_pure_bdipole(tm, bunch, bdict[1])
  end
end

@inline function thin_pure_bquadrupole(tm::Linear, bunch, bm2)
  K1L = get_thin_strength(bm2, 0, bunch.Brho_ref)
  mx, my = LinearTracking.linear_thin_quad_matrices(K1L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, 0, nothing, nothing))
end

# Ignore higher order multipoles with quadrupole: 
@inline thin_bquadrupole(tm::Linear, bunch, bdict) = thin_pure_bquadrupole(tm, bunch, bdict[2])

# Ignore higher order multipoles: do nothing for thin case (0 length)
@inline thin_pure_bmultipole(tm::Linear, bunch, bmn) = KernelCall()
@inline thin_bmultipole(tm::Linear, bunch, bdict) = KernelCall()


# === Thick elements === #
@inline function drift(tm::Linear, bunch, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  return KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2))
end

@inline function thick_pure_bsolenoid(tm::Linear, bunch, bm0, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  Ks = get_thick_strength(bm0, L, bunch.Brho_ref)
  mxy = LinearTracking.linear_solenoid_matrix(Ks, L)
  return KernelCall(LinearTracking.linear_coast!, (mxy, L/gamma_0^2, nothing, nothing))
end

@inline function thick_pure_bdipole(tm::Linear, bunch, bm1, L)
  # Sophia: this a thick corrector coil
  # In Fortran Bmad is g == 0 but dg != 0. 
  # In SciBmad this is g == 0 but K0 != 0.
  error("Undefined for tracking method $tm")
end

@inline function thick_bdipole(tm::Linear, bunch, bdict, L)
  if haskey(bdict, 2)
    # Sophia: this is a thick corrector coil with a quadrupole term
    # In Fortran Bmad is g == 0 but dg != 0 and K1 != 0
    # In SciBmad this is g == 0 but K0 != 0 and K1 != 0
    gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
    K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
    K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref) 
    mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(K0, L, gamma_0; K1 = K1)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
  else # ignore higher order multipoles
    return thick_pure_bdipole(tm, bunch, bdict[1], L)
  end
end

@inline function thick_pure_bquadrupole(tm::Linear, bunch, bm2, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  mx, my = LinearTracking.linear_quad_matrices(K1, L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, L/gamma_0^2, nothing, nothing))
end

# Ignore higher order multipoles with quadrupole:
@inline thick_bquadrupole(tm::Linear, bunch, bdict, L) = thick_pure_bquadrupole(tm, bunch, bdict[2], L)

# Ignore higher order multipoles: drift
@inline thick_pure_bmultipole(tm::Linear, bunch, bmn, L) = drift(tm, bunch, L)
@inline thick_bmultipole(tm::Linear, bunch, bdict, L) = drift(tm, bunch, L)

# =========== BENDING ELEMENTS ============= #
# Again, "bend" means coordinate system bends, nothing wrt the field.
# SciBmad likely will not support "thin" bends

@inline function thick_bend_no_field(tm::Linear, bunch, bendparams, L)
  # Sophia: this has NO FIELD!
  # In Fortran Bmad it is like setting dg == -g
  error("Undefined for tracking method $tm")
end

@inline function thick_bend_pure_bdipole(tm::Linear, bunch, bendparams, bm1, L)     
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(K0, L, gamma_0; g = bendparams.g, e1 = bendparams.e1, e2 = bendparams.e2)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
end

@inline function thick_bend_bdipole(tm::Linear, bunch, bendparams, bdict, L)   
  if haskey(bdict, 2)
    # Sophia: this is a thick combined function magnet
    # In Fortran Bmad is g != 0, dg != 0, K1 != 0
    # In SciBmad this is g != 0, K0 != 0, K1 != 0
    gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
    K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
    K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref)
    mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(K0, L, gamma_0; g = bendparams.g, K1 = K1, e1=bendparams.e1, e2=bendparams.e2)
    return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
  else # ignore higher order multipoles
    return thick_bend_pure_bdipole(tm, bunch, bendparams, bdict[1], L)
  end
end

@inline function thick_bend_pure_bquadrupole(tm::Linear, bunch, bendparams, bm2, L) 
  # Sophia: this is a quadrupole with a g, I think your code should be able to handle this
  # In Fortran Bmad it would be like dg == -g, K1 != 0.
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(0, L, gamma_0; g = bendparams.g, K1 = K1, e1=bendparams.e1, e2=bendparams.e2)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
end

# Ignore higher order multipoles:
@inline thick_bend_bquadrupole(tm::Linear, bunch, bendparams, bdict, L) = thick_bend_pure_bquadrupole(tm, bunch, bendparams, bdict[2], L)
# Ignore higher order multipoles: treat like bend no field
@inline thick_bend_pure_bmultipole(tm::Linear, bunch, bendparams, bmn, L) = thick_bend_no_field(tm, bunch, bendparams, L)
@inline thick_bend_bmultipole(tm::Linear, bunch, bendparams, bdict, L) = thick_bend_no_field(tm, bunch, bendparams, L)
