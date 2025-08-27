_track!(
  i,
  coords::Coords,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm::SciBmadStandard;
  kwargs...
) = _track!(i, coords, bunch, ele, SplitIntegration(); kwargs...)