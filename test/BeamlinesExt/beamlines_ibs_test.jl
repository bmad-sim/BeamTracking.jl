using Random 

@testset "IBS" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  drift = Drift(L = 2.0, 
  tracking_method = Symplectic(ibs_damping_on = true, ibs_fluctuations_on = true))
  line = Beamline([drift], species_ref = Species("electron"), E_ref = 18e9)

  v0 = randn(10,6)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line; use_explicit_SIMD = !(VERSION < v"1.11" && Sys.ARCH == :x86_64))

  # The kick on live particles does not depend on a lost particle's coordinates
  line = Beamline([Drift(L = 1.0, tracking_method = Symplectic(ibs_damping_on = true))],
                  species_ref = Species("proton"), E_ref = 1e9)
  sig = [1e-4 1e-5 1e-4 1e-5 1e-3 1e-4]
  v0 = randn(50, 6) .* sig
  b_live = Bunch(v = copy(v0), weight = 1e4, species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  b_lost = Bunch(v = [v0; 20 .* sig], state = [fill(BeamTracking.STATE_ALIVE, 50); BeamTracking.STATE_LOST],
                 weight = 1e4, species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b_live, line; use_explicit_SIMD = !(VERSION < v"1.11" && Sys.ARCH == :x86_64))
  track!(b_lost, line; use_explicit_SIMD = !(VERSION < v"1.11" && Sys.ARCH == :x86_64))
  @test b_lost.coords.v[1:50, [2, 4, 6]] - v0[:, [2, 4, 6]] ≈ b_live.coords.v[:, [2, 4, 6]] - v0[:, [2, 4, 6]]

  # Six live particles cannot define a nonsingular 6D covariance
  b_few = Bunch(v = copy(v0), state = [fill(BeamTracking.STATE_ALIVE, 6); fill(BeamTracking.STATE_LOST, 44)],
                weight = 1e4, species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  @test_throws ErrorException track!(b_few, line; use_explicit_SIMD = !(VERSION < v"1.11" && Sys.ARCH == :x86_64))
end
