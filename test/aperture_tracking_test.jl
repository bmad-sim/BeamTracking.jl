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

v0 = [
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
  @test bunch.coords.state == [STATE_LOST, STATE_ALIVE, STATE_LOST, STATE_LOST, STATE_LOST, STATE_LOST]

  bunch = Bunch(copy(v2))
  BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_aperture_elliptical!, (1.0, 2.0, 3.0, 5.0)))
  @test bunch.coords.state == [STATE_LOST, STATE_ALIVE]

  test_matrix(v0, KernelCall(BeamTracking.track_aperture_rectangular!, (1.0, 2.0, 3.0, 5.0)))
  test_matrix(v0, KernelCall(BeamTracking.track_aperture_elliptical!, (1.0, 2.0, 3.0, 5.0)))
end