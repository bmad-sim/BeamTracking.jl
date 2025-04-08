module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA
using BeamTracking: soaview, get_N_particle, MAX_TEMPS, calc_gamma, launch!
import BeamTracking: track!, make_track_chain

# Define my own custom tracking method
# This should be in Beamlines itself so default to this?
struct SciBmadStandard end



include("linear.jl")


end