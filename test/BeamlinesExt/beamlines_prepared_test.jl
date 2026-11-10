@testset "Prepared tracking" begin
  e = Species("electron")
  R = Beamlines.pc_to_R(e, 18e9)
  v0 = [1e-3 1e-4 -1e-3 2e-4 1e-3 1e-4; -2e-3 0 1e-3 -1e-4 0 -2e-4; 1e-4 0 1e-4 0 0 0]
  q0 = repeat([1.0 0 0 0], 3)
  # RF (the reference time grows from turn to turn), slices, a time dependent strength, a
  # misalignment, an aperture, a thin multipole and a marker
  ring(p_over_q_ref) = Beamline([
    LatElement(L=1.0, Kn1=0.3, tracking_method=MatrixKick(order=4), x_offset=1e-4,
               x1_limit=-0.1, x2_limit=0.1, y1_limit=-0.1, y2_limit=0.1),
    LatElement(L=2.0, n_slices=2, n_steps=2, Kn1=-0.25, tracking_method=DriftKick(order=4)),
    LatElement(L=0.5, voltage=1e6, rf_frequency=5e8, zero_phase=PhaseRef.AboveTransition,
               tracking_method=SolenoidKick(order=4)),
    LatElement(L=0.5, Kn1=0.3*(1 + 1e5*Time()), tracking_method=DriftKick(order=2), n_steps=2),
    LatElement(L=0.0, Kn2L=0.5),
    LatElement(L=0.0),
  ]; p_over_q_ref=p_over_q_ref, species_ref=e)

  # Tracks bitwise like track!(bunch, bl) over many turns
  function test_prepared(bl, mk; turns=20, kw...)
    b1 = mk(); b2 = mk()
    pbl = btbl.prepare(b2, bl; notify=false, kw...)
    for t in 1:turns
      track!(b1, bl; notify=false, kw...)
      track!(b2, pbl)
    end
    @test b1.coords.v == b2.coords.v
    @test isequal(b1.coords.q, b2.coords.q)
    @test b1.coords.state == b2.coords.state
    @test b1.t_ref == b2.t_ref
    @test b1.p_over_q_ref == b2.p_over_q_ref
    return pbl
  end
  bl = ring(R)
  pbl = test_prepared(bl, () -> Bunch(copy(v0); p_over_q_ref=R, species=e))
  @test length(pbl.chains) == 7
  test_prepared(bl, () -> Bunch(Float32.(v0); p_over_q_ref=R, species=e))
  test_prepared(bl, () -> Bunch(copy(v0), copy(q0); p_over_q_ref=R, species=e))
  test_prepared(bl, () -> Bunch(copy(v0); p_over_q_ref=R, species=e); rf_on=false)
  # Another species than the reference (z is re-referenced to the reference species' time)
  test_prepared(bl, () -> Bunch(copy(v0); p_over_q_ref=R, species=Species("muon")))
  # SaganCavity (without a reference energy change, so that it can be tracked repeatedly)
  bl6 = lat2.beamlines[4]
  R6 = bl6.p_over_q_ref
  test_prepared(bl6, () -> Bunch(copy(v0), copy(q0); species=species, p_over_q_ref=R6))

  # Only Bunches like the one prepared with; not supported yet
  @test_throws ErrorException track!(Bunch(Float32.(v0); p_over_q_ref=R, species=e), pbl)
  @test_throws ErrorException track!(Bunch(copy(v0); p_over_q_ref=2R, species=e), pbl)
  @test_throws ErrorException btbl.prepare(Bunch(copy(v0)), ring(R*(1 + 1e3*Time())); notify=false)
end
