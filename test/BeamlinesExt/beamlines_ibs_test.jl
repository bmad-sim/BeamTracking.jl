using Random 

@testset "IBS" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  drift = Drift(L = 2.0, 
  tracking_method = Yoshida(ibs_damping_on = true, ibs_fluctuations_on = true))
  line = Beamline([drift], species_ref = Species("electron"), E_ref = 18e9)

  v0 = randn(10,6)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)
end
