module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions
using Beamlines: isactive, BitsLineElement
using BeamTracking: soaview, get_N_particle, calc_gamma, calc_p0c, runkernels!,
                    @makekernel, BunchView, KernelCall, KernelChain, push
import BeamTracking: track!, C_LIGHT, chargeof, massof


include("utils.jl")

function track!(
  bunch::Bunch, 
  ele::LineElement; 
  kwargs...
)
  b = BunchView(bunch)
  @noinline _track!(nothing, b, bunch, ele, ele.tracking_method; kwargs...)
  return bunch
end

# Indicies array instead of nothing
# this could be used to actually mask alive particles, specify indicies
# Would also allow you to do mix of outer and inner loop too, doing a sub-bunch of 
# particles in parallel

@makekernel fastgtpsa=false function outer_track!(i, b::BunchView, bunch::Bunch, bl::Beamline)
  for j in 1:length(bl.line)
    @inbounds ele = bl.line[j]
    @noinline _track!(i, b, bunch, ele, ele.tracking_method)
  end
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  outer_particle_loop::Bool=false,
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end

  check_Brho(bl.Brho_ref, bunch)

  if !outer_particle_loop
    for ele in bl.line
      track!(bunch, ele; kwargs...)
    end
  else
    kc = (KernelCall(outer_track!, (bunch, bl)),)
    launch!(BunchView(bunch), kc; kwargs...)
  end

  return bunch
end


function track!(
  bunch::Bunch, 
  bbl::BitsBeamline{TM}; 
  outer_particle_loop::Bool=false
) where {TM}

  if length(bbl.params) == 0
    return bunch
  end
  
  check_Brho(NaN, bunch)

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
            _track!(nothing, BunchView(bunch), bunch, ele, TM)
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
        _track!(nothing, BunchView(bunch), bunch, ele, TM)
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
  outer_particle_loop::Bool=false
) where {TM<:Beamlines.MultipleTrackingMethods}
  error("BitsBeamline tracking including different tracking methods per element not implemented yet")
end
#=
function _track!(
  i,
  b::BunchView,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm::Any;
  kwargs...
)
  error("Tracking method $tm is not defined!")
end
=#
include("unpack.jl")
include("linear.jl")
include("exact.jl")


end