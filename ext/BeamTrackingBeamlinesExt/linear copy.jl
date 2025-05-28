
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]


# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  v,
  work,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::Linear;
  kwargs...
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L

  # Function barrier
  linear_universal!(i, v, work, bunch, L, bm, bp, ma; kwargs...)
end

@inline function get_thick_strength(bm, L, Brho_ref)
  s = bm.strength
  if !bm.normalized
    s /= Brho_ref
  end
  if bm.integrated
    if L == 0
      error("LineElement length is zero; cannot computed non-integrated strength")
    end
    s /= L
  end
  return s
end

@inline function get_thin_strength(bm, L, Brho_ref)
  s = bm.strength
  if !bm.normalized
    s /= Brho_ref
  end
  if !bm.integrated
    s *= L
  end
  return s
end

# Step 2: Push particles through -----------------------------------------
function linear_universal!(
  i, 
  v, 
  work,
  bunch,
  L, 
  bmultipoleparams, 
  bendparams, 
  alignmentparams;
  kwargs...
  #=
  groupsize::Union{Nothing,Integer}=nothing, #backend isa CPU ? floor(Int,REGISTER_SIZE/sizeof(eltype(v))) : 256 
  multithread_threshold::Integer=0,
  use_KA::Bool=!(get_backend(v) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=false
  =#
) 
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  if !isactive(bmultipoleparams) # Drift
    if isactive(bendparams)
      error("Linear tracking requires BendParams.g == BMultipoleParams.K0")
    end
    runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2; kwargs...)
  elseif haskey(bmultipoleparams.bdict, 0) # Solenoid
    if any(t -> t >= 1, keys(bmultipoleparams.bdict))
      error("Linear tracking does not support combined solenoid + other multipole magnets")
    end
    if isactive(bendparams)
      error("Linear tracking does not currently support solenoid with bending")
    end
    if L == 0
      error("Thin solenoid not supported yet")
    end

    Ks = get_thick_strength(bmultipoleparams.bdict[0], L, bunch.Brho_ref)

    mxy = LinearTracking.linear_solenoid_matrix(Ks, L)
    runkernel!(LinearTracking.linear_coast!, i, v, work, mxy, L/gamma_0^2, nothing, nothing; kwargs...)
  elseif haskey(bmultipoleparams.bdict, 1) # Bend
    if !isactive(bendparams)
      error("Linear tracking requires BendParams.g ≈ BMultipoleParams.K0")
    end
    if L == 0
      error("Thin bend not supported yet")
    end
    if any(t -> t == 0 || t > 1, keys(bmultipoleparams.bdict))
      error("Combined function bend tracking not implemented yet")
    end

    K0 = get_thick_strength(bmultipoleparams.bdict[1], L, bunch.Brho_ref)

    if !(K0 ≈ bendparams.g)
      error("Linear tracking requires BendParams.g ≈ BMultipoleParams.K0")
    end
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, gamma_0, bendparams.e1, bendparams.e2)
    runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56, d, t; kwargs...)
  elseif haskey(bmultipoleparams.bdict, 2) # Quadrupole
    if isactive(bendparams)
      error("For Linear combined function magnet tracking, both the K0 multipole and BendParams must be set")
    end
    if L == 0
      K1L = get_thin_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
      mx, my = LinearTracking.linear_thin_quad_matrices(K1L)
    else
      K1 = get_thick_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
      mx, my = LinearTracking.linear_quad_matrices(K1, L)
    end
    runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, L/gamma_0^2, nothing, nothing; kwargs...)
  else # Drift for higher-order multipoles
    runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2; kwargs...)
  end

  return v
end
