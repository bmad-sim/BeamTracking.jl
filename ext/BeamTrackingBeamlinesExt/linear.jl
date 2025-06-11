
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

  if isactive(bendparams) #geometric bend
    if haskey(bmultipoleparams.bdict, 0) 
      error("Linear tracking does not currently support solenoid with bending")
    end 
    #bend
    if haskey(bmultipoleparams.bdict,1)
      if L == 0
        error("Thin bend not supported yet")
      end
      if any(t -> t == 0 || t > 2, keys(bmultipoleparams.bdict)) 
        error("Linear tracking does not support bend tracking including any other multipole except a quadrupole")
      end
      K0 = get_thick_strength(bmultipoleparams.bdict[1], L, bunch.Brho_ref)
      if !haskey(bmultipoleparams.bdict, 2)
        K1 = nothing
      else
        K1 = get_thick_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref) 
      end
      mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(K0, L, gamma_0; g=bendparams.g, K1=K1, e1=bendparams.e1, e2=bendparams.e2)
      runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56, d, t; kwargs...)
    else
      error("Geometric bend specified without field K0")
    end
  else
    if !isactive(bmultipoleparams) #drift
      runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2; kwargs...)
    elseif haskey(bmultipoleparams.bdict, 0) # Solenoid
      if any(t -> t >= 1, keys(bmultipoleparams.bdict))
        error("Linear tracking does not support combined solenoid + other multipole magnets")
      end
      if L == 0
        error("Thin solenoid not supported yet")
      end
      Ks = get_thick_strength(bmultipoleparams.bdict[0], L, bunch.Brho_ref)
      mxy = LinearTracking.linear_solenoid_matrix(Ks, L)
      runkernel!(LinearTracking.linear_coast!, i, v, work, mxy, L/gamma_0^2; kwargs...)
    elseif haskey(bmultipoleparams.bdict, 1) # Bend
      if L == 0
        error("Thin bend not supported yet")
      end
      if any(t -> t == 0 || t > 2, keys(bmultipoleparams.bdict)) 
        error("Linear tracking does not support bend tracking including any other multipole except a quadrupole")
      end
      K0 = get_thick_strength(bmultipoleparams.bdict[1], L, bunch.Brho_ref)
      if !haskey(bmultipoleparams.bdict, 2)
        K1 = nothing
      else
        K1 = get_thick_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref) 
      end
      mx, my, r56, d, t = LinearTracking.linear_dipole_matrices(K0, L, gamma_0; g=nothing, K1=K1, e1=bendparams.e1, e2=bendparams.e2)
      runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56, d, t; kwargs...)
    elseif haskey(bmultipoleparams.bdict, 2) # Quadrupole
      if L == 0
        K1L = get_thin_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
        mx, my = LinearTracking.linear_thin_quad_matrices(K1L)
      else
        K1 = get_thick_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
        mx, my = LinearTracking.linear_quad_matrices(K1, L)
      end
      runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, L/gamma_0^2; kwargs...)
    else # Drift for higher-order multipoles
      runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2; kwargs...)
    end 
  end
  return v
end
