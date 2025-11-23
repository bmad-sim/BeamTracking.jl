_track!(
  i,
  coords::Coords,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm::SciBmadStandard,
  ramp_without_rf;
  kwargs...
) = _track!(i, coords, bunch, ele, SplitIntegration(), ramp_without_rf; kwargs...)