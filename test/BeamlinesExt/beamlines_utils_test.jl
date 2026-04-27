const btbl = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)

@testset "BeamlinesUtils" begin
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1), Species("electron")) ≈ 0.1 + pi
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.Accelerating), Species("positron")) ≈ 0.1
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.BelowTransition), Species("positron")) ≈ 0.1 + 0.5*pi
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.AboveTransition), Species("electron")) ≈ 0.1 - 0.5*pi + pi
end