using BeamTracking, Test

v1 = [ 
       0.0  0.0  0.1  0.0  0.0  0.0
     ]

tilde_m = 0.5
beta_0 = 1 / sqrt(1 + tilde_m^2)
gamsqr_0 = 1 + 1 / tilde_m^2

# Signature is:
#   @makekernel fastgtpsa=true function track_alignment_straight!(i, coords::Coords, beta_0, gamsqr_0, tilde_m,
#      entering, x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)


@testset "AlignmentKernel" begin
  bunch = Bunch(copy(v1))
  BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_alignment_straight!, (beta_0, gamsqr_0, tilde_m,
                          true, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, +1, 2.0)))
  println(v1)
  println(bunch.coords.v)
  @test bunch.coords.v == v1
end

