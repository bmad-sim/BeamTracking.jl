include("../lattices/sagan_cavity_lat.jl")
import BeamTracking as BT

v1 = [0.0  0.0  0.0  0.0  0.0  0.0]
elec = Species("electron")
R0 = BT.E_to_R(elec, 1e6)

#@testset "SaganCavity" begin
  b1 = Bunch(deepcopy(v1), species=elec, p_over_q_ref = R0)
  track!(b1, sc1)
  println(b1.coords.v)
#end



# Test spin
# Test multipoles
# Test zero length
# Test phi0
# standing/traveling waves
# dE_ref
# differing n_steps
# differing n_cell
# non-zero starting orbit
