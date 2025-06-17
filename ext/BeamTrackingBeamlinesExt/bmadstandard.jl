@inline function bend_fringe(tm::BmadStandard, bunch, bendparams, bm1, L, upstream)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  ge = upstream ? K0*tan(bendparams.e1) : K0*tan(bendparams.e2)
  return KernelCall(BmadStandardTracking.linear_dipole_hard_edge_fringe!, (ge,))
end

@inline function thick_bend_pure_bdipole(tm::BmadStandard, bunch, bendparams, bm1, L)
  g = bendparams.g
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  me2 = bendparams.e2 ≈ 0 ? 0 : K0*tan(bendparams.e2)
  return KernelCall(BmadStandardTracking.combined_func!, (L, g, K0, 0, tilde_m))
end

@inline function thick_bend_bdipole(tm::BmadStandard, bunch, bendparams, bdict, L)
  g = bendparams.g
  K0 = get_thick_strength(bdict[1], L, bunch.Brho_ref)
  K1 = get_thick_strength(bdict[2], L, bunch.Brho_ref)
  tilde_m = massof(bunch.species)/calc_p0c(bunch.species, bunch.Brho_ref)
  me2 = bendparams.e2 ≈ 0 ? 0 : K0*tan(bendparams.e2)
  return KernelCall(BmadStandardTracking.combined_func!, (L, g, K0, K1, tilde_m))
end