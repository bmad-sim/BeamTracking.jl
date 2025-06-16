# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  b::BunchView,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm;
  kwargs...
)
  # Unpack the line element (type unstable)
  L = ele.L
  ap = ele.AlignmentParams
  bp = ele.BendParams
  bm = ele.BMultipoleParams
  # bc BitsLineElement does not support PatchParams yet
  pp = ele isa BitsLineElement ? nothing : ele.PatchParams

  # Function barrier
  universal!(i, b, tm, bunch, L, ap, bp, bm, pp; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  i, 
  b,
  tm,
  bunch,
  L, 
  alignmentparams, 
  bendparams,
  bmultipoleparams,
  patchparams;
  kwargs...
) 
  kc = KernelChain(Val{1}())
  if isactive(alignmentparams)
    if isactive(patchparams)
      error("Tracking through a LineElement containing both PatchParams and AlignmentParams is undefined")
    end
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, true)))
  end

  if isactive(patchparams) 
    # Patch
    if isactive(bendparams) || isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both PatchParams and BendParams/BMultipoleParams not currently defined")
    else
      kc = push(kc, @inline(pure_patch(tm, bunch, patchparams, L)))
    end
  elseif isactive(bendparams)
    # Bend
    if !isactive(bmultipoleparams) 
      # Bend no field
      kc = push(kc, @inline(bend_no_field(tm, bunch, bendparams, L)))
    else
      bdict = bmultipoleparams.bdict
      n_multipoles = length(keys(bdict))
      if haskey(bdict, 0) # Bend-solenoid
        if n_multipoles == 1
          # Pure bend-solenoid
          kc = push(kc, @inline(bend_pure_bsolenoid(tm, bunch, bendparams, bdict[0], L)))
        else
          # Bend-solenoid with other multipoles of order > 0
          kc = push(kc, @inline(bend_bsolenoid(tm, bunch, bendparams, bdict, L)))
        end
      elseif haskey(bdict, 1) # Bend-dipole
        if n_multipoles == 1
          # Pure bend-dipole
          kc = push(kc, @inline(bend_pure_bdipole(tm, bunch, bendparams, bdict[1], L)))
        else
          # Bend-dipole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bdipole(tm, bunch, bendparams, bdict, L)))
        end
      elseif haskey(bdict, 2) # Bend-quadrupole
        if n_multipoles == 1
          # Pure bend-quadrupole
          kc = push(kc, @inline(bend_pure_bquadrupole(tm, bunch, bendparams, bdict[2], L)))
        else
          # Bend-quadrupole with other multipoles of order > 1
          kc = push(kc, @inline(bend_bquadrupole(tm, bunch, bendparams, bdict, L)))
        end
      else # Bend-multipole
        if n_multipoles == 1
          # Pure bend-multipole
          kc = push(kc, @inline(bend_pure_bmultipole(tm, bunch, bendparams, first(values(bdict)), L)))
        else
          # Bend-multipole with other multipoles of order > 2
          kc = push(kc, @inline(bend_bmultipole(tm, bunch, bendparams, bdict, L)))
        end
      end
    end
  elseif isactive(bmultipoleparams)
    # BMultipole
    bdict = bmultipoleparams.bdict
    n_multipoles = length(keys(bdict))
    if haskey(bdict, 0) # Solenoid
      if n_multipoles == 1
        # Pure solenoid
        kc = push(kc, @inline(pure_bsolenoid(tm, bunch, bdict[0], L)))
      else
        # Solenoid with other multipoles of order > 0
        kc = push(kc, @inline(bsolenoid(tm, bunch, bdict, L)))
      end
    elseif haskey(bdict, 1) # Dipole without bend
      if n_multipoles == 1
        # Pure dipole
        kc = push(kc, @inline(pure_bdipole(tm, bunch, bdict[1], L)))
      else
        # Dipole with other multipoles of order > 1
        kc = push(kc, @inline(bdipole(tm, bunch, bdict, L)))
      end
    elseif haskey(bdict, 2) # Quadrupole
      if n_multipoles == 1
        # Pure quadrupole
        kc = push(kc, @inline(pure_bquadrupole(tm, bunch, bdict[2], L)))
      else
        # Quadrupole with other multipoles of order > 1
        kc = push(kc, @inline(bquadrupole(tm, bunch, bdict, L)))
      end
    else # Higher order multipole
      if n_multipoles == 1
        # Pure multipole
        kc = push(kc, @inline(pure_bmultipole(tm, bunch, first(values(bdict)), L)))
      else
        # Multipole with other multipoles of order > 2
        kc = push(kc, @inline(bmultipole(tm, bunch, bdict, L)))
      end
    end
  elseif !(L ≈ 0)
    kc = push(kc, @inline(drift(tm, bunch, L)))
  end
    
  if isactive(alignmentparams)
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, false)))
  end

  runkernels!(i, b, kc; kwargs...)
  return nothing
end

# === Coordinate system transformations === #
@inline misalign(tm, bunch, alignmentparams, in) = error("Undefined for tracking method $tm")
@inline pure_patch(tm, bunch, patchparams, L)    = error("Undefined for tracking method $tm")

# === Straight Elements === #
# "Pure" means only ONE SINGLE multipole.
# When "pure" is not present, it means that at least one HIGHER ORDER
# multipole exists.
@inline thin_pure_bsolenoid(tm, bunch, bm0)   = error("Undefined for tracking method $tm")
@inline thin_bsolenoid(tm, bunch, bdict)      = error("Undefined for tracking method $tm")
@inline thin_pure_bdipole(tm, bunch, bm1)     = error("Undefined for tracking method $tm")
@inline thin_bdipole(tm, bunch, bdict)        = error("Undefined for tracking method $tm")
@inline thin_pure_bquadrupole(tm, bunch, bm2) = error("Undefined for tracking method $tm")
@inline thin_bquadrupole(tm, bunch, bdict)    = error("Undefined for tracking method $tm")
@inline thin_pure_bmultipole(tm, bunch, bmn)  = error("Undefined for tracking method $tm")
@inline thin_bmultipole(tm, bunch, bdict)     = error("Undefined for tracking method $tm")

@inline thick_pure_bsolenoid(tm, bunch, bm0, L)   = error("Undefined for tracking method $tm")
@inline thick_bsolenoid(tm, bunch, bdict, L)      = error("Undefined for tracking method $tm")
@inline thick_pure_bdipole(tm, bunch, bm1, L)     = error("Undefined for tracking method $tm")
@inline thick_bdipole(tm, bunch, bdict, L)        = error("Undefined for tracking method $tm")
@inline thick_pure_bquadrupole(tm, bunch, bm2, L) = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, bunch, bdict, L)    = error("Undefined for tracking method $tm")
@inline thick_pure_bmultipole(tm, bunch, bmn, L)  = error("Undefined for tracking method $tm")
@inline thick_bmultipole(tm, bunch, bdict, L)     = error("Undefined for tracking method $tm")


# === Elements with curving coordinate system "bend" === #
# "Bend" means ONLY a coordinate system curvature through the element.
# It does NOT IMPLY ANY DIPOLE FIELD! Bend specifies if the integration 
# path is curving but does not IMPACT THE PHYSICS inside.
# "Pure" means only ONE SINGLE MULTIPOLE
# When "pure" is not present, it means that at least one higher order 
# multipole exists.

# SciBmad will probably not support thin bends ever but I leave them here for now
@inline thin_bend_no_field(tm, bunch, bendparams)              = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bsolenoid(tm, bunch, bendparams, bm0)   = error("Undefined for tracking method $tm")
@inline thin_bend_bsolenoid(tm, bunch, bendparams, bdict)      = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bdipole(tm, bunch, bendparams, bm1)     = error("Undefined for tracking method $tm")
@inline thin_bend_bdipole(tm, bunch, bendparams, bdict)        = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bquadrupole(tm, bunch, bendparams, bm2) = error("Undefined for tracking method $tm")
@inline thin_bend_bquadrupole(tm, bunch, bendparams, bdict)    = error("Undefined for tracking method $tm")
@inline thin_bend_pure_bmultipole(tm, bunch, bendparams, bmn)  = error("Undefined for tracking method $tm")
@inline thin_bend_bmultipole(tm, bunch, bendparams, bdict)     = error("Undefined for tracking method $tm")

@inline thick_bend_no_field(tm, bunch, bendparams, L)              = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)   = error("Undefined for tracking method $tm")
@inline thick_bend_bsolenoid(tm, bunch, bendparams, bdict, L)      = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bdipole(tm, bunch, bendparams, bm1, L)     = error("Undefined for tracking method $tm")
@inline thick_bend_bdipole(tm, bunch, bendparams, bdict, L)        = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L) = error("Undefined for tracking method $tm")
@inline thick_bend_bquadrupole(tm, bunch, bendparams, bdict, L)    = error("Undefined for tracking method $tm")
@inline thick_bend_pure_bmultipole(tm, bunch, bendparams, bmn, L)  = error("Undefined for tracking method $tm")
@inline thick_bend_bmultipole(tm, bunch, bendparams, bdict, L)     = error("Undefined for tracking method $tm")


# === Elements thin vs thick check === #
@inline pure_bsolenoid(tm, bunch, bm0, L)   = L ≈ 0 ? thin_pure_bsolenoid(tm, bunch, bm0)   : thick_pure_bsolenoid(tm, bunch, bm0, L)      
@inline bsolenoid(tm, bunch, bdict, L)      = L ≈ 0 ? thin_bsolenoid(tm, bunch, bdict)      : thick_bsolenoid(tm, bunch, bdict, L)         
@inline pure_bdipole(tm, bunch, bm1, L)     = L ≈ 0 ? thin_pure_bdipole(tm, bunch, bm1)     : thick_pure_bdipole(tm, bunch, bm1, L)          
@inline bdipole(tm, bunch, bdict, L)        = L ≈ 0 ? thin_bdipole(tm, bunch, bdict)        : thick_bdipole(tm, bunch, bdict, L)             
@inline pure_bquadrupole(tm, bunch, bm2, L) = L ≈ 0 ? thin_pure_bquadrupole(tm, bunch, bm2) : thick_pure_bquadrupole(tm, bunch, bm2, L)        
@inline bquadrupole(tm, bunch, bdict, L)    = L ≈ 0 ? thin_bquadrupole(tm, bunch, bdict)    : thick_bquadrupole(tm, bunch, bdict, L)           
@inline pure_bmultipole(tm, bunch, bmn, L)  = L ≈ 0 ? thin_pure_bmultipole(tm, bunch, bmn)  : thick_pure_bmultipole(tm, bunch, bmn, L)                   
@inline bmultipole(tm, bunch, bdict, L)     = L ≈ 0 ? thin_bmultipole(tm, bunch, bdict)     : thick_bmultipole(tm, bunch, bdict, L)                      
@inline bend_no_field(tm, bunch, bendparams, L)              = L ≈ 0 ? thin_bend_no_field(tm, bunch, bendparams)              : thick_bend_no_field(tm, bunch, bendparams, L)
@inline bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)   = L ≈ 0 ? thin_bend_pure_bsolenoid(tm, bunch, bendparams, bm0)   : thick_bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)      
@inline bend_bsolenoid(tm, bunch, bendparams, bdict, L)      = L ≈ 0 ? thin_bend_bsolenoid(tm, bunch, bendparams, bdict)      : thick_bend_bsolenoid(tm, bunch, bendparams, bdict, L)         
@inline bend_pure_bdipole(tm, bunch, bendparams, bm1, L)     = L ≈ 0 ? thin_bend_pure_bdipole(tm, bunch, bendparams, bm1)     : thick_bend_pure_bdipole(tm, bunch, bendparams, bm1, L)          
@inline bend_bdipole(tm, bunch, bendparams, bdict, L)        = L ≈ 0 ? thin_bend_bdipole(tm, bunch, bendparams, bdict)        : thick_bend_bdipole(tm, bunch, bendparams, bdict, L)             
@inline bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L) = L ≈ 0 ? thin_bend_pure_bquadrupole(tm, bunch, bendparams, bm2) : thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)        
@inline bend_bquadrupole(tm, bunch, bendparams, bdict, L)    = L ≈ 0 ? thin_bend_bquadrupole(tm, bunch, bendparams, bdict)    : thick_bend_bquadrupole(tm, bunch, bendparams, bdict, L)           
@inline bend_pure_bmultipole(tm, bunch, bendparams, bmn, L)  = L ≈ 0 ? thin_bend_pure_bmultipole(tm, bunch, bendparams, bmn)  : thick_bend_pure_bmultipole(tm, bunch, bendparams, bmn, L)                   
@inline bend_bmultipole(tm, bunch, bendparams, bdict, L)     = L ≈ 0 ? thin_bend_bmultipole(tm, bunch, bendparams, bdict)     : thick_bend_bmultipole(tm, bunch, bendparams, bdict, L)                      

