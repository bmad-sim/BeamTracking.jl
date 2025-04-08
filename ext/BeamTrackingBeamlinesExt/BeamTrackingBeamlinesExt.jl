module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA
using BeamTracking: soaview, get_N_particle, MAX_TEMPS, calc_gamma, launch!
import BeamTracking: track!

# Define my own custom tracking method
# This should be in Beamlines itself so default to this?
struct SciBmadStandard end

function track!(bunch::Bunch, bl::Beamline, bunch0::Bunch; bunch_work::Bunch=deepcopy(bunch0), work=nothing)
  soaview(bunch_work) .= soaview(bunch0)
  bunch_work.species = bunch0.species
  bunch_work.Brho_0 = bunch0.Brho_0

  # Preallocate the temporaries if not provided
  if isnothing(work)
    n_temps = 0
    for ele in bl.line
      cur_n_temps = MAX_TEMPS(ele.tracking_method) 
      n_temps = cur_n_temps > n_temps ? cur_n_temps : n_temps
    end
    N_particle = get_N_particle(bunch) 
    work = zeros(eltype(bunch.v), N_particle, n_temps)
  end

  b1 = bunch_work
  b2 = bunch
  for ele in bl.line
    track!(b2, ele, b1; work=work)
    b1, b2 = b2, b1 # pointer swap
  end

  # if odd number of elements then bunch contains result
  # if even then bunch_work contains the result, so one final step to copy:
  if iseven(length(bl.line))
    bunch.species = bunch_work.species
    bunch.Brho_0 = bunch_work.Brho_0
    v = soaview(bunch)
    v_work = soaview(bunch_work)
    # use this instead of copy because then I handle mutability of TPS numbers safely
    foreach(i->@FastGTPSA!(v[i] = v_work[i]), 1:get_N_particle(bunch))
  end
  return bunch
end


function track!(bunch::Bunch, ele::LineElement, bunch0::Bunch; work=nothing)
  # Dispatch on the tracking method:
  return track!(bunch, ele, bunch0, ele.tracking_method; work=work)
end

include("linear.jl")


end