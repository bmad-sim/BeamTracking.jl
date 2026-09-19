@testset "BeamlinesDoNotUse" begin
  species = Species("electron")
  v0 = [1e-3 2e-4 -5e-4 1e-4 1e-3 1e-3]

  function track_ele(ele)
    b = Bunch(copy(v0), species=species)
    track!(b, Beamline([ele], species_ref=species, E_ref=1e9))
    return b
  end
  same(b1, b2) = b1.coords.v ≈ b2.coords.v && b1.coords.state == b2.coords.state

  # AlignmentParams
  q_ref = Quadrupole(L=0.5, Kn1=0.3)
  q = Quadrupole(L=0.5, Kn1=0.3, x_offset=1e-3, tilt=0.1)
  @test !same(track_ele(q), track_ele(q_ref))
  q.do_not_use = [:AlignmentParams]
  @test same(track_ele(q), track_ele(q_ref))

  # BMultipoleParams: quadrupole becomes a drift
  q = Quadrupole(L=0.5, Kn1=0.3, do_not_use=[:BMultipoleParams])
  @test same(track_ele(q), track_ele(Drift(L=0.5)))

  # ApertureParams
  q = Quadrupole(L=0.5, Kn1=0.3, x1_limit=-1e-4, x2_limit=1e-4, y1_limit=-1e-4, y2_limit=1e-4)
  @test track_ele(q).coords.state[1] != STATE_ALIVE
  q.do_not_use = [:ApertureParams]
  @test same(track_ele(q), track_ele(q_ref))

  # BendParams: bend becomes a straight element with the same multipoles, and
  # alignment is done as for a straight element
  sb = SBend(L=2.0, g=0.01, Kn0=0.01, x_offset=1e-3, y_rot=1e-3)
  straight = LineElement(L=2.0, Kn0=0.01, x_offset=1e-3, y_rot=1e-3)
  @test !same(track_ele(sb), track_ele(straight))
  sb.do_not_use = [:BendParams]
  @test same(track_ele(sb), track_ele(straight))
  sb.do_not_use = [:BendParams, :BMultipoleParams, :AlignmentParams]
  @test same(track_ele(sb), track_ele(Drift(L=2.0)))

  # RFParams: cavity becomes a drift
  cav = RFCavity(L=1.0, voltage=1e6, rf_frequency=5e8)
  @test !same(track_ele(cav), track_ele(Drift(L=1.0)))
  cav.do_not_use = [:RFParams]
  @test same(track_ele(cav), track_ele(Drift(L=1.0)))

  # PatchParams
  pa = Patch(dx=1e-3)
  @test !same(track_ele(pa), track_ele(Marker()))
  pa.do_not_use = [:PatchParams]
  @test same(track_ele(pa), track_ele(Marker()))

  # All instances of an element in a Beamline use the do_not_use list of the element
  q = Quadrupole(L=0.5, Kn1=0.3, x_offset=1e-3)
  bl = Beamline([q, Drift(L=1.0), q], species_ref=species, E_ref=1e9)
  bl_ref = Beamline([q_ref, Drift(L=1.0), q_ref], species_ref=species, E_ref=1e9)
  b1 = Bunch(copy(v0), species=species); track!(b1, bl)
  b_ref = Bunch(copy(v0), species=species); track!(b_ref, bl_ref)
  @test !same(b1, b_ref)
  q.do_not_use = [:AlignmentParams]
  b2 = Bunch(copy(v0), species=species); track!(b2, bl)
  @test same(b2, b_ref)

  # Invalid symbols added directly to the list are caught when tracking
  push!(q.do_not_use, :AlignmentParam)
  @test_throws ErrorException track!(Bunch(copy(v0), species=species), bl)
end
