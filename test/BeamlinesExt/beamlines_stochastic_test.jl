using Random 

@testset "Stochastic radiation" begin
  Random.seed!(0)

  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  bend = SBend(g = 0.01, L = 2.0, 
  tracking_method = Yoshida(order = 2, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true))
  line = Beamline([bend], species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)

  v0 = [0.01 0.02 0.03 0.04 0.05 0.06]
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.048892166423554095 
  0.02117516850989616 
  0.10556186816301837 
  0.040000228166203126 
  0.0476104945848721 
  0.06000604532567636]'

  bend.tracking_method = Yoshida(order = 2, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.0488918794021746 
  0.021174170337409773 
  0.1055618678196179 
  0.03999901918537246 
  0.047610501975114294 
  0.05997400853912311]'

  bend.tracking_method = Yoshida(order = 4, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889191059360618 
  0.02117425155266759 
  0.10556186829371136 
  0.039999117112352 
  0.04761050112770837 
  0.05997659826202305]'

  bend.tracking_method = Yoshida(order = 4, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889188935716087 
  0.021174111144523146 
  0.10556186753391111 
  0.03999898748624612 
  0.04761050156666369 
  0.05997317490679352]'

  bend.tracking_method = Yoshida(order = 6, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889204991364092 
  0.021174378272207466 
  0.10556186803046733 
  0.0399990778319358 
  0.047610497477820424 
  0.0599755625943892]'

  bend.tracking_method = Yoshida(order = 6, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)
  
  @test b0.coords.v ≈
  [0.04889192415056395 
  0.021174272439899 
  0.10556186809803797 
  0.0399991405912617 
  0.047610500761089414 
  0.05997722303511241]'

  bend.tracking_method = Yoshida(order = 8, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈
  [0.048892222474638355 
  0.02117491594142489 
  0.10556186814738264 
  0.039999748582531414 
  0.047610492940657556 
  0.059993338012956456]'

  bend.tracking_method = Yoshida(order = 8, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈
  [0.0488920515629509 
  0.02117473434154834 
  0.10556186792526312
  0.03999967905433932 
  0.04761049754544777 
  0.059991498329934695]'
end
 
