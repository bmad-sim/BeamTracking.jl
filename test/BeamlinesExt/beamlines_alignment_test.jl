include("../lattices/alignment_lat.jl")

v1 = [ 
       0.0  0.0  1.0  0.0  0.0  0.0
     ]

@testset "Alignment" begin
  b1 = Bunch(copy(v1), species=Species("electron"), R_ref = 1.0)
  track!(b1, bline_d1)
  b2 = Bunch(copy(v1), species=Species("electron"), R_ref = 1.0)
  track!(b2, bline_d2)
  ## @test b1.v ≈ b2.v
end

