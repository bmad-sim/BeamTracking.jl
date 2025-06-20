using CUDA, BeamTracking, Beamlines, PhysicalConstants

import  PhysicalConstants.CODATA2022: c_0 as c, m_e as m

include("../../test/lattices/esr.jl") # Beamline symbol is "ring"
# Currently only Linear tracking is supported, enable it for each element
foreach(t -> t.tracking_method = Linear(), ring.line)

N_particle = 100
b0 = Bunch(N_particle)

bitsring = BitsBeamline(ring)
track!(b0, bitsring)