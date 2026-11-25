@testset "BeamlinesStraightBendKick" begin
  # A dipole with multipoles in a straight reference frame (g_ref = 0), e.g. a rectangular bend
  # placed between patches: the multipoles are Maxwellian, so no "curved reference system" warning
  species = Species("electron")
  v0 = [1e-3 2e-4 -5e-4 1e-4 1e-3 1e-3]
  ele = SBend(L=2.0, g_ref=0.0, Kn0=0.01, Kn1=1e-3, Kn2=0.1, Ks3=0.5, tracking_method=BendKick(order=6), n_steps=4)
  bl = Beamline([ele], species_ref=species, E_ref=18e9)
  b = Bunch(copy(v0), species=species)
  @test_logs track!(b, bl)
  @test all(isfinite, b.coords.v)
  @test b.coords.state[1] == BeamTracking.STATE_ALIVE

  # The multipoles act in the straight frame: with them off, the same element is the pure dipole
  b1 = Bunch(copy(v0), species=species)
  track!(b1, Beamline([SBend(L=2.0, g_ref=0.0, Kn0=0.01, tracking_method=BendKick())], species_ref=species, E_ref=18e9))
  @test !(b1.coords.v ≈ b.coords.v)
end
