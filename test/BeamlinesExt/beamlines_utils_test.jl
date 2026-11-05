const btbl = Base.get_extension(BeamTracking, :BeamTrackingBeamlinesExt)

@testset "BeamlinesUtils" begin
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1), Species("electron")) ≈ 0.1 + pi
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.Accelerating), Species("positron")) ≈ 0.1
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.BelowTransition), Species("positron")) ≈ 0.1 + 0.5*pi
  @test btbl.rf_phi0_calc(RFParams(phi0 = 0.1, zero_phase = PhaseRef.AboveTransition), Species("electron")) ≈ 0.1 - 0.5*pi + pi
end

@testset "BeamlinesReference" begin
  e = Species("electron")
  p = Species("positron")
  R = BeamTracking.E_to_R(e, 1e9) # Signed: negative for electrons
  v0 = [1e-3 2e-4 -5e-4 1e-4 1e-3 1e-3]

  # Species and energy unset in the Bunch are inherited from the Beamline
  b = Bunch(copy(v0))
  ref = btbl.resolve_reference!(b, e, R, false)
  @test b.species == e && b.p_over_q_ref ≈ R
  @test ref.species_ref == e && ref.p_over_q_ref ≈ R && ref.field_scale == 1

  # ... and unset in the Beamline are inherited from the Bunch
  b = Bunch(copy(v0); species=e, p_over_q_ref=1.1*R)
  ref = btbl.resolve_reference!(b, nothing, nothing, false)
  @test ref.species_ref == e && ref.p_over_q_ref == 1.1*R && b.p_over_q_ref == 1.1*R

  # Both set but different: the Bunch keeps its energy here, and is ramped when tracking
  b = Bunch(copy(v0); species=e, p_over_q_ref=1.1*R)
  ref = btbl.resolve_reference!(b, e, R, false)
  @test ref.p_over_q_ref == R && b.p_over_q_ref == 1.1*R

  # Errors
  @test_throws ErrorException btbl.resolve_reference!(Bunch(copy(v0)), nothing, R, false)                         # No species
  @test_throws ErrorException track!(Bunch(copy(v0); species=e), Beamline([Drift(L=1.0)]); notify=false)          # No energy
  @test_throws ErrorException btbl.resolve_reference!(Bunch(copy(v0); species=e, p_over_q_ref=-R), e, R, false)   # Wrong sign
  @test_throws ErrorException btbl.resolve_reference!(Bunch(copy(v0); species=p, p_over_q_ref=-R), e, nothing, false) # Fields have no physical value

  # Opposite charge: same magnetic rigidity magnitude, normalized strengths flip sign
  b = Bunch(copy(v0); species=p)
  ref = btbl.resolve_reference!(b, e, R, false)
  @test ref.field_scale == -1 && ref.p_over_q_ref == -R && b.p_over_q_ref == -R
  kn, ks = btbl.get_strengths(LatElement(L=2.0, Kn1=0.3).BMultipoleParams, 2.0, -R, -1)
  @test kn[1] ≈ -0.3
  kn, ks = btbl.get_strengths(LatElement(L=2.0, Bn1=0.3).BMultipoleParams, 2.0, -R, -1)
  @test kn[1] ≈ 0.3/(-R)

  # A positron through a Beamline designed for electrons (fields held fixed) is the same as
  # a positron through a Beamline designed for positrons with the same fields
  bl_e = Beamline([Quadrupole(L=0.5, Kn1=0.3), Drift(L=1.0), SBend(L=1.0, g=0.01, Kn0=0.01)]; species_ref=e, p_over_q_ref=R)
  bl_p = Beamline([Quadrupole(L=0.5, Kn1=-0.3), Drift(L=1.0), SBend(L=1.0, g=0.01, Kn0=-0.01)]; species_ref=p, p_over_q_ref=-R)
  b1 = Bunch(copy(v0); species=p); track!(b1, bl_e; notify=false)
  b2 = Bunch(copy(v0); species=p); track!(b2, bl_p; notify=false)
  @test b1.coords.v ≈ b2.coords.v
  @test b1.t_ref ≈ b2.t_ref

  # Different mass: a proton through a Beamline designed for deuterons at the same rigidity sees
  # the same fields, but z and t_ref are referenced to the deuteron (design) clock
  d = Species("deuteron")
  pr = Species("proton")
  R_d = BeamTracking.E_to_R(d, 3e9)
  bl_d = Beamline([Quadrupole(L=0.5, Kn1=0.3), Drift(L=1.0)]; species_ref=d, p_over_q_ref=R_d)
  bl_pr = Beamline([Quadrupole(L=0.5, Kn1=0.3), Drift(L=1.0)]; species_ref=pr, p_over_q_ref=R_d)
  b1 = Bunch(copy(v0); species=pr); track!(b1, bl_d; notify=false)
  b2 = Bunch(copy(v0); species=pr); track!(b2, bl_pr; notify=false)
  @test b1.coords.v[:, [1,2,3,4,6]] ≈ b2.coords.v[:, [1,2,3,4,6]]
  t_d = 1.5/BeamTracking.R_to_v(d, R_d)
  t_pr = 1.5/BeamTracking.R_to_v(pr, R_d)
  @test b1.t_ref ≈ t_d && b2.t_ref ≈ t_pr
  beta = (1 + v0[6])/sqrt((1 + v0[6])^2 + 1/BeamTracking.R_to_beta_gamma(pr, R_d)^2)
  @test b1.coords.v[1,5] ≈ b2.coords.v[1,5] + beta*C_LIGHT*(t_d - t_pr)
end
