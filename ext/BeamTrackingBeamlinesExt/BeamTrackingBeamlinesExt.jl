module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays
using Beamlines: isactive, BitsLineElement
using BeamTracking: soaview, get_N_particle, calc_gamma, launch!, runkernel!
import BeamTracking: track!, MAX_TEMPS

# Specify a MAX_TEMPS for SciBmadStandard
MAX_TEMPS(::SciBmadStandard) = 1

include("utils.jl")

function track!(
  bunch::Bunch, 
  ele::LineElement; 
  work=zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(ele.tracking_method)),
)
  return _track!(nothing, soaview(bunch), work, bunch, ele, ele.tracking_method)
end

# Indicies array instead of nothing
# this could be used to actually mask alive particles, specify indicies
# Would also allow you to do mix of outer and inner loop too, doing a sub-bunch of 
# particles in parallel

@inline function outer_track!(i, v, work, bunch, bl::Beamline)
  for j in 1:length(bl.line)
    @inbounds ele = bl.line[j]
    @noinline _track!(i, v, work, bunch, ele, ele.tracking_method)
  end
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  work=get_work(bunch, bl), 
  outer_particle_loop::Bool=false
)
  if length(bl.line) == 0
    return bunch
  end

  check_Brho(bl.Brho_ref, bunch)

  if !outer_particle_loop
    for ele in bl.line
      track!(bunch, ele; work=work)
    end
  else
    launch!(outer_track!, soaview(bunch), work, bunch, bl)
  end

  return bunch
end


function track!(
  bunch::Bunch, 
  bbl::BitsBeamline{TM}; 
  work=get_work(bunch, bbl), 
  outer_particle_loop::Bool=false
) where {TM}

  if length(bbl.params) == 0
    return bunch
  end

  if !outer_particle_loop
    if !isnothing(bbl.rep)
      i = 1 
      while i <= length(bbl.params)
        repeat_count = bbl.rep[i]
        start_i = i
        for _ in 1:repeat_count
          i = start_i
          while true
            ele = BitsLineElement(bbl, i)
            _track!(nothing, soaview(bunch), work, bunch, ele, TM)
            i += 1
            if i > length(bbl.rep) || bbl.rep[i] != 0
              break
            end
          end
        end
      end
    else
      for i in 1:length(bbl.params)
        ele = BitsLineElement(bbl, i)
        _track!(nothing, soaview(bunch), work, bunch, ele, TM)
      end
    end
  else
    error("outer_particle_loop tracking for BitsBeamline not implemented yet")
  end

  return bunch
end

function track!(
  bunch::Bunch, 
  bbl::BitsBeamline{TM}; 
  work=get_work(bunch, bbl), 
  outer_particle_loop::Bool=false
) where {TM<:Beamlines.MultipleTrackingMethods}
  error("BitsBeamline tracking including different tracking methods per element not implemented yet")
end



include("linear.jl")


end