using Random 

@testset "IBS" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  drift = Drift(L = 2.0, 
  tracking_method = Yoshida(ibs_damping_on = true, ibs_fluctuations_on = true))
  line = Beamline([drift], species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)

  v0 = randn(10,6)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors, and equivalent with IBS on vs. off
  drift.tracking_method = Yoshida(ibs_damping_on = true, ibs_fluctuations_on = true)
  b0 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors
  drift.tracking_method = Yoshida(ibs_damping_on = false, ibs_fluctuations_on = false)
  b02 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b02, line)

  @test norm(normTPS.(b0.coords.v - b02.coords.v)) ≈ 0
end