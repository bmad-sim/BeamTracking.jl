using Random 

@testset "Stochastic radiation" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  bend = SBend(g = 0.01, L = 2.0, 
  tracking_method = Yoshida(order = 2, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true))
  line = Beamline([bend], species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)

  v0 = [0.01 0.02 0.03 0.04 0.05 0.06]
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 2, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 4, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 4, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 6, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 6, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 8, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Yoshida(order = 8, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Now just check SIMD , if it doesn't bug out
  bend.tracking_method = Yoshida(order = 8, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(rand(10,6), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors, and equivalent with fluctuations on vs. off
  bend.tracking_method = Yoshida(order = 2, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors
  bend.tracking_method = Yoshida(order = 2, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = false)
  b02 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b02, line)

  @test norm(normTPS.(b0.coords.v - b02.coords.v)) ≈ 0
end