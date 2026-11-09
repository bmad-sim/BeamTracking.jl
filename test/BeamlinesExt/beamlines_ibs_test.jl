using Random 

@testset "IBS" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  drift = Drift(L = 2.0, 
  tracking_method = DriftKick(), ibs_model = IBSModel.Gaussian)
  line = Beamline([drift], species_ref = Species("electron"), E_ref = 18e9)

  v0 = randn(10,6)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Lost particles (frozen far away) do not count: the live particles see the same IBS as
  # without them, also with weights (damping only, so deterministic)
  drift.ibs_model = IBSModel.GaussianDampingOnly
  line = Beamline([drift], species_ref = Species("electron"), E_ref = 100e6)
  v1 = randn(100,6) .* [1e-4 1e-5 1e-5 1e-6 1e-3 1e-4]
  v_lost = fill(1.0, 5, 6)
  state = vcat(fill(BeamTracking.STATE_ALIVE, 100), fill(BeamTracking.STATE_LOST, 5))
  for (w, w_lost) in ((nothing, nothing), (1e8, 1e8), (fill(1e8, 100), fill(5e9, 5)))
    b1 = Bunch(v = copy(v1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref, weight = w)
    b2 = Bunch(v = vcat(v1, v_lost), state = copy(state), species = line.species_ref, p_over_q_ref = line.p_over_q_ref,
      weight = w isa Vector ? vcat(w, w_lost) : w)
    track!(b1, line)
    track!(b2, line)
    @test b1.coords.v[:,6] != v1[:,6] # IBS changed the momenta
    @test b2.coords.v[1:100,:] ≈ b1.coords.v rtol=1e-12
    @test b2.coords.v[101:end,:] == v_lost
  end
end
