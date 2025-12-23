const btbl = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)

@testset "BeamlinesUtils" begin
  @test btbl.rf_phi0(RFParams(phi0 = 0.1)) ≈ 0.1
  @test btbl.rf_phi0(RFParams(phi0 = 0.1, zero_phase = PhaseReference.Accelerating)) ≈ 0.1
  @test btbl.rf_phi0(RFParams(phi0 = 0.1, zero_phase = PhaseReference.BelowTransition)) ≈ 0.1 + 0.5*pi
  @test btbl.rf_phi0(RFParams(phi0 = 0.1, zero_phase = PhaseReference.AboveTransition)) ≈ 0.1 - 0.5*pi
end