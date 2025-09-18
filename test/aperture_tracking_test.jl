v1 = [ 
       0.0  0.0  6.0  0.0  0.0  0.0
       1.9  0.8  4.0  0.0  0.0  0.0
       2.1 -0.8  4.0  0.0  0.0  0.0
       0.0  0.0  4.0  0.0  0.0  0.0
       0.0  0.0  2.0  0.0  0.0  0.0
       1.5  0.0 -6.0  0.0  0.0  0.0
     ]

v2 = [ 
       1.9  0.0  4.9  0.0  0.0  0.0
       1.8  0.9  4.2  0.0  0.0  0.0
     ]

m_unit = [
        1.0 0.0 0.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0 0.0 0.0
        0.0 0.0 1.0 0.0 0.0 0.0
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
     ]

@testset "ApertureKernel" begin
  bunch = Bunch(copy(v1))
  BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_aperture_rectangular!, (1.0, 2.0, 3.0, 5.0)))
  @test bunch.coords.state == [STATE_LOST_POS_Y, STATE_ALIVE, STATE_LOST_POS_X, STATE_LOST_NEG_X, STATE_LOST_NEG_Y, STATE_LOST_NEG_Y]

  bunch = Bunch(copy(v2))
  BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_aperture_elliptical!, (1.0, 2.0, 3.0, 5.0)))
  @test bunch.coords.state == [STATE_LOST_POS_Y, STATE_ALIVE]

  test_matrix(m_unit, KernelCall(BeamTracking.track_aperture_rectangular!, (1.0, 2.0, 3.0, 5.0)))
  test_matrix(m_unit, KernelCall(BeamTracking.track_aperture_elliptical!, (1.0, 2.0, 3.0, 5.0)))
end