module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays, KernelAbstractions
using Beamlines: isactive, deval, unsafe_getparams, o2i, BitsBeamline, BitsLineElement, isnullspecies
using BeamTracking: get_N_particle, R_to_beta_gamma, R_to_gamma, R_to_pc, R_to_v, beta_gamma_to_v, runkernels!,
                    @makekernel, Coords, KernelCall, KernelChain, push, TimeDependentParam, RefState, launch!
import BeamTracking: track!, C_LIGHT, chargeof, massof


include("utils.jl")

function track!(
  bunch::Bunch, 
  ele::LineElement; 
  t_ref::Ref=Ref{eltype(bunch.coords.v)}(0),
  kwargs...
)
  coords = bunch.coords
  @noinline _track!(nothing, coords, bunch, t_ref, ele, ele.tracking_method; kwargs...)
  return bunch
end

# Indicies array instead of nothing
# this could be used to actually mask alive particles, specify indicies
# Would also allow you to do mix of outer and inner loop too, doing a sub-bunch of 
# particles in parallel

@makekernel fastgtpsa=false function outer_track!(i, b::Coords, bunch::Bunch, bl::Beamline)
  for j in 1:length(bl.line)
    @inbounds ele = bl.line[j]
    @noinline _track!(i, b, bunch, ele, ele.tracking_method)
  end
end

function track!(
  bunch::Bunch, 
  bl::Beamline; 
  t_ref::Ref=Ref{eltype(bunch.coords.v)}(0),
  outer_particle_loop::Bool=false,
  kwargs...
)
  if length(bl.line) == 0
    return bunch
  end

  check_bl_bunch!(bl, bunch)

  if !outer_particle_loop
    for ele in bl.line
      track!(bunch, ele; t_ref=t_ref, kwargs...)
    end
  else
    kc = (KernelCall(outer_track!, (bunch, bl)),)
    launch!(bunch.coords, kc; kwargs...)
  end

  return bunch
end


function track!(
  bunch::Bunch, 
  bbl::BitsBeamline{TM}; 
  t_ref::Ref=Ref{eltype(bunch.coords.v)}(0),
  outer_particle_loop::Bool=false
) where {TM}

  if length(bbl.params) == 0
    return bunch
  end
  
  check_R_ref!(nothing, bunch)

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
            _track!(nothing, bunch.coords, bunch, t_ref, ele, TM)
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
        _track!(nothing, bunch.coords, bunch, t_ref, ele, TM)
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
  b::Coords,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm::Any;
  kwargs...
)
  error("Tracking method $tm is not defined!")
end

=#

include("aperture.jl")
include("alignment.jl")
include("unpack.jl")
include("scibmadstandard.jl")
include("linear.jl")
include("exact.jl")
include("integration.jl")

end