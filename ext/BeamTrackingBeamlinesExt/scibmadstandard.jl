_track!(
  i,
  b::BunchView,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm::SciBmadStandard;
  kwargs...
) = _track!(i, b, bunch, ele, SplitIntegration(); kwargs...)