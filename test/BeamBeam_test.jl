#(i, coords::Coords, p0c, charge,
	# beta_0, tilde_m,
	# sig_x_strong, sig_y_strong, sig_z_strong, N_particles,
	# n_slices, z_offset)
# using BeamTracking,
#     Beamlines,
#     JET,
#     BenchmarkTools,
#     GTPSA,
#     StaticArrays,
#     ReferenceFrameRotations,
#     SIMD,
#     KernelAbstractions
# import KernelAbstractions: @kernel
# using BeamTracking: Coords, KernelCall, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST
# using Beamlines: isactive
#include("../src/kernels/BeamBeam_tracking.jl")
xi  = [ 1.000]
pxi = [ 0.100]
yi  = [ 0.000]
pyi = [ 0.000]
zi  = [ -1.000]
pzi = [ 0.900]

xifinal = [1.0078]
pxifinal = [0.102958]
yifinal = [0.0000]
pyifinal = [0.0000]
zifinal = [-1.00004]
pzifinal = [0.900077]

@testset "BeamBeamTracking" begin
    @testset "BeamBeam" begin
        v = [ xi pxi yi pyi zi pzi ]
        bunch = Bunch(v)
        BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (100000000,1e11,1,0.1,0.1,0,1e14,1,0)))
        print(bunch.coords)
        @test v[:,BeamTracking.XI]  ≈  xifinal (rtol=5.e-4)
        @test v[:,BeamTracking.YI]  ≈  yifinal (rtol=5.e-4)
        @test v[:,BeamTracking.ZI]  ≈  zifinal (rtol=5.e-4)
        @test v[:,BeamTracking.PXI] == pxifinal (rtol=5.e-4)
        @test v[:,BeamTracking.PYI] == pyifinal (rtol=5.e-4)
        @test v[:,BeamTracking.PZI] == pzifinal (rtol=5.e-4)
    end
end


