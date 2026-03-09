using Test
using AtomicAndPhysicalConstants
import BeamTracking as BT

@testset "Miscellaneous" begin
  @test BT.R_to_v(Species("He--"), -0.3) ≈ 1.445522537423045e7
end