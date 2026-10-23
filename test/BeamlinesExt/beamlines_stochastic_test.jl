using Random 

@testset "Stochastic radiation" begin
  Random.seed!(0)

  # Here just check that they don't bug out
  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  bend = SBend(g = 0.01, L = 2.0, 
  tracking_method = Symplectic(order = 2, n_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true))
  line = Beamline([bend], species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)

  v0 = [0.01 0.02 0.03 0.04 0.05 0.06]
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 2, n_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 4, n_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 4, n_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 6, n_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 6, n_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 8, n_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  bend.tracking_method = Symplectic(order = 8, n_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Now just check SIMD , if it doesn't bug out
  bend.tracking_method = Symplectic(order = 8, n_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(rand(10,6), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors, and equivalent with fluctuations on vs. off
  bend.tracking_method = Symplectic(order = 2, n_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  # Track TPSA make sure no errors
  bend.tracking_method = Symplectic(order = 2, n_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = false)
  b02 = Bunch(vars(D1), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b02, line)

  @test norm(normTPS.(b0.coords.v - b02.coords.v)) ≈ 0

  # A tilted bend acts as rotating into the bend frame, tracking the untilted bend and rotating
  # back. With the same random numbers this holds for the radiation fluctuations as well.
  function dipole(x, y, s, t, p)
    g, Kn0 = p
    potential = (0.0, 0.0, 0.0, -Kn0*(x + g*x^2/2))
    jac = (0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           -Kn0*(1 + g*x), 0.0, 0.0, 0.0)
    return potential, jac
  end
  function rot(v, a)
    c = cos(a)
    s = sin(a)
    return [c*v[1]+s*v[3] c*v[2]+s*v[4] -s*v[1]+c*v[3] -s*v[2]+c*v[4] v[5] v[6]]
  end
  function track_seeded(kw, v)
    tm = Symplectic(order = 2, n_steps = 2, radiation_damping_on = true, radiation_fluctuations_on = true)
    line = Beamline([LineElement(; L = 2.0, tracking_method = tm, kw...)],
                    species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)
    b = Bunch(copy(v), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
    Random.seed!(1)
    track!(b, line)
    return b.coords.v
  end
  tilt = 0.4
  for kw in ((g = 0.01,), (g = 0.01, Kn1 = 0.1),
             (g_ref = 0.01, four_potential = dipole, four_potential_params = (0.01, 0.01), four_potential_normalized = true))
    @test track_seeded((; tilt_ref = tilt, kw...), v0) ≈ rot(track_seeded(kw, rot(v0, tilt)), -tilt) atol=1e-14
  end
  # A separator whose field components (En0, Es0) rotate its field by the same angle
  E = 1e7
  @test track_seeded((En0 = E*cos(tilt), Es0 = -E*sin(tilt)), v0) ≈ rot(track_seeded((En0 = E,), rot(v0, tilt)), -tilt) atol=1e-14
end
