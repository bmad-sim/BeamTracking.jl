include("../lattices/alignment_lat.jl")

v1 = [0.1  0.2  0.3  0.4  0.5  0.6]

@testset "Alignment" begin
  b1 = Bunch(copy(v1), species=Species("electron"), R_ref = 1.0)
  track!(b1, bline_d1)
  b2 = Bunch(copy(v1), species=Species("electron"), R_ref = 1.0)
  track!(b2, bline_d2)

  # Test that misalignments do not effect tracking through a drift.
  @test b1.coords.v ≈ b2.coords.v
end

