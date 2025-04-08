function track!(
  bunch::Bunch, 
  ele::LineElement, 
  ::Linear; 
  work=nothing
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L
  Brho_ref = ele.Brho_ref

  # Function barrier (this function is now fully compiled)
  # For some reason, inlining this is faster/zero allocs
  # copy-paste is slower and so is @noinline so I guess LLVM is 
  # doing some kind of constant propagation while inlining this?
  return @inline _track_linear!(bunch, ma, bm, bp, L, Brho_ref; work=work)
end

function _track_linear!(
  bunch::Bunch, 
  ma::Union{AlignmentParams,Nothing},
  bm::Union{BMultipoleParams,Nothing}, 
  bp::Union{BendParams,Nothing},
  L, 
  Brho_ref;
  work=nothing # =zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(Linear())),
) 
  v = soaview(bunch)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_0)

  if !isnothing(bp)
    error("Bend tracking not implemented yet")
  end

  if !isnothing(ma)
    #chain = merge(chain, (ExactTracking.misalign!, ma.x_offset, ma.y_offset, -1))
  end

  if isnothing(bm) || length(bm.bdict) == 0 # Drift
    launch!(LinearTracking.linear_drift!, v, nothing, L, L/gamma_0^2)
    #chain = merge(chain, (LinearTracking.linear_drift!, L, gamma_0))
  else
    if length(bm.bdict) > 1 || !haskey(bm.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end

    bm1 = bm.bdict[2]
    if bm1.tilt != 0
      error("Currently multipole tilts not supported")
    end

    # Tracking should work with B-fields. The lattice has no understanding of a particle beam
    # or the different types of particles shooting through

    # The Bunch struct itself can store the normalization factor of its own phase space coordinates 
    # but that may in general be different from the lattice.

    # So we always just get the B-fields from the lattice
    B1 = bm1.strength
    if bm1.normalized
      B1 *= Brho_ref
    end
    if bm1.integrated
      B1 /= L
    end

    # Now get K1 from the bunch itself:
    K1 = B1/bunch.Brho_0

    mx, my = LinearTracking.linear_quad_matrices(K1, L)
    r56 = L/gamma_0^2 
    wq = isnothing(work) ? zeros(eltype(bunch.v), get_N_particle(bunch), 1) : work
    launch!(LinearTracking.linear_coast_uncoupled!, v, wq, mx, my, r56)
    #chain = merge(chain, (LinearTracking.linear_coast_uncoupled!, mx, my, r56))
  end

  if !isnothing(ma)
    #chain = merge(chain, (ExactTracking.misalign!, ma.x_offset, ma.y_offset, -1))
  end

  #launch!(chain, v, work)
  return bunch
end