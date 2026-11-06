using Test
using AtomicAndPhysicalConstants
import BeamTracking as BT

@testset "Miscellaneous" begin
  @test Beamlines.R_to_v(Species("He--"), -0.3) ≈ 1.4442639106654184e7

  # Energy conversions must be consistent for any charge state (|Z| != 1 included).
  @testset "Energy conversions $name" for (name, R) in (("electron", -1.0), ("proton", 1.0),
                                    ("He++", 1.0), ("He--", -0.3), ("C+6", 1.0), ("C+6", 20.0))
    sp = Species(name)
    E = Beamlines.R_to_E(sp, R)
    beta_gamma = Beamlines.R_to_beta_gamma(sp, R)
    @test Beamlines.R_to_gamma(sp, R) ≈ E/massof(sp)
    @test Beamlines.R_to_gamma(sp, R) ≈ sqrt(1 + beta_gamma^2)
    @test Beamlines.R_to_v(sp, R) ≈ BT.beta_gamma_to_v(beta_gamma)
    @test Beamlines.R_to_v(sp, R) ≈ Beamlines.E_to_v(sp, E)
    @test Beamlines.R_to_v(sp, R) < BT.C_LIGHT
    @test Beamlines.E_to_R(sp, E) ≈ R
    @test Beamlines.pc_to_R(sp, Beamlines.R_to_pc(sp, R)) ≈ R
  end
end
