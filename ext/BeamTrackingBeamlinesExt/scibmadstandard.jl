_track!(
  i,
  coords::Coords,
  bunch::Bunch,
  t_ref::Ref,
  ele::Union{LineElement,BitsLineElement}, 
  tm::SciBmadStandard;
  kwargs...
) = _track!(i, coords, bunch, t_ref, ele, SplitIntegration(); kwargs...)