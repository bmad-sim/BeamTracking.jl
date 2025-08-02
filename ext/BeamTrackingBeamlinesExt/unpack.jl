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
  L = ele.L # Automatically calls deval (element-level get)
  ap = deval(ele.AlignmentParams)
  bp = deval(ele.BendParams)
  bm = deval(ele.BMultipoleParams)
  pp = deval(ele.PatchParams)
  dp = deval(ele.ApertureParams)

  # Function barrier
  universal!(i, b, tm, bunch, L, ap, bp, bm, pp, dp; kwargs...)
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
  patchparams,
  apertureparams;
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
@inline thin_pure_bsolenoid(tm, bunch, bm0)           = error("Undefined for tracking method $tm")
@inline thin_bsolenoid(tm, bunch, bmultipoleparams)   = error("Undefined for tracking method $tm")
@inline thin_pure_bdipole(tm, bunch, bm1)             = error("Undefined for tracking method $tm")
@inline thin_bdipole(tm, bunch, bmultipoleparams)     = error("Undefined for tracking method $tm")
@inline thin_pure_bquadrupole(tm, bunch, bm2)         = error("Undefined for tracking method $tm")
@inline thin_bquadrupole(tm, bunch, bmultipoleparams) = error("Undefined for tracking method $tm")
@inline thin_pure_bmultipole(tm, bunch, bmk)          = error("Undefined for tracking method $tm")
@inline thin_bmultipole(tm, bunch, bmultipoleparams)  = error("Undefined for tracking method $tm")

@inline thick_pure_bsolenoid(tm, bunch, bm0, L)           = error("Undefined for tracking method $tm")
@inline thick_bsolenoid(tm, bunch, bmultipoleparams, L)   = error("Undefined for tracking method $tm")
@inline thick_pure_bdipole(tm, bunch, bm1, L)             = error("Undefined for tracking method $tm")
@inline thick_bdipole(tm, bunch, bmultipoleparams, L)     = error("Undefined for tracking method $tm")
@inline thick_pure_bquadrupole(tm, bunch, bm2, L)         = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, bunch, bmultipoleparams, L) = error("Undefined for tracking method $tm")
@inline thick_pure_bmultipole(tm, bunch, bmk, L)          = error("Undefined for tracking method $tm")
@inline thick_bmultipole(tm, bunch, bmultipoleparams, L)  = error("Undefined for tracking method $tm")


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
@inline pure_bsolenoid(tm, bunch, bm0, L)                            = L ≈ 0 ? thin_pure_bsolenoid(tm, bunch, bm0)                            : thick_pure_bsolenoid(tm, bunch, bm0, L)      
@inline bsolenoid(tm, bunch, bmultipoleparams, L)                    = L ≈ 0 ? thin_bsolenoid(tm, bunch, bmultipoleparams)                    : thick_bsolenoid(tm, bunch, bmultipoleparams, L)         
@inline pure_bdipole(tm, bunch, bm1, L)                              = L ≈ 0 ? thin_pure_bdipole(tm, bunch, bm1)                              : thick_pure_bdipole(tm, bunch, bm1, L)          
@inline bdipole(tm, bunch, bmultipoleparams, L)                      = L ≈ 0 ? thin_bdipole(tm, bunch, bmultipoleparams)                      : thick_bdipole(tm, bunch, bmultipoleparams, L)             
@inline pure_bquadrupole(tm, bunch, bm2, L)                          = L ≈ 0 ? thin_pure_bquadrupole(tm, bunch, bm2)                          : thick_pure_bquadrupole(tm, bunch, bm2, L)        
@inline bquadrupole(tm, bunch, bmultipoleparams, L)                  = L ≈ 0 ? thin_bquadrupole(tm, bunch, bmultipoleparams)                  : thick_bquadrupole(tm, bunch, bmultipoleparams, L)           
@inline pure_bmultipole(tm, bunch, bmk, L)                           = L ≈ 0 ? thin_pure_bmultipole(tm, bunch, bmk)                           : thick_pure_bmultipole(tm, bunch, bmk, L)                   
@inline bmultipole(tm, bunch, bmultipoleparams, L)                   = L ≈ 0 ? thin_bmultipole(tm, bunch, bmultipoleparams)                   : thick_bmultipole(tm, bunch, bmultipoleparams, L)                      
@inline bend_no_field(tm, bunch, bendparams, L)                      = L ≈ 0 ? thin_bend_no_field(tm, bunch, bendparams)                      : thick_bend_no_field(tm, bunch, bendparams, L)
@inline bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)           = L ≈ 0 ? thin_bend_pure_bsolenoid(tm, bunch, bendparams, bm0)           : thick_bend_pure_bsolenoid(tm, bunch, bendparams, bm0, L)      
@inline bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)   = L ≈ 0 ? thin_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams)   : thick_bend_bsolenoid(tm, bunch, bendparams, bmultipoleparams, L)         
@inline bend_pure_bdipole(tm, bunch, bendparams, bm1, L)             = L ≈ 0 ? thin_bend_pure_bdipole(tm, bunch, bendparams, bm1)             : thick_bend_pure_bdipole(tm, bunch, bendparams, bm1, L)          
@inline bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)     = L ≈ 0 ? thin_bend_bdipole(tm, bunch, bendparams, bmultipoleparams)     : thick_bend_bdipole(tm, bunch, bendparams, bmultipoleparams, L)             
@inline bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)         = L ≈ 0 ? thin_bend_pure_bquadrupole(tm, bunch, bendparams, bm2)         : thick_bend_pure_bquadrupole(tm, bunch, bendparams, bm2, L)        
@inline bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L) = L ≈ 0 ? thin_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams) : thick_bend_bquadrupole(tm, bunch, bendparams, bmultipoleparams, L)           
@inline bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)          = L ≈ 0 ? thin_bend_pure_bmultipole(tm, bunch, bendparams, bmk)          : thick_bend_pure_bmultipole(tm, bunch, bendparams, bmk, L)                   
@inline bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)  = L ≈ 0 ? thin_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams)  : thick_bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)                      

