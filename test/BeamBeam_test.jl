using GTPSA
using FiniteDifferences
d_z = Descriptor(6, 1)
b0 = Bunch(collect(transpose(@vars(d_z))), R_ref=-0.0017045090263411496, species=Species("electron"))

xi  = [ 1.000]
pxi = [ 0.100]
yi  = [ 2.000]
pyi = [ 0.200]
zi  = [ -0.1]
pzi = [ 0.900]

xifinal = [9.9701282840875e-01]
pxifinal = [1.0011665425195e-01]
yifinal = [1.9940256568175e00]
pyifinal = [2.0023330850389e-01]
zifinal = [-9.9213481467006e-02]
pzifinal = [9.0009482033428e-01]

xi2  = [ 1.000, 2.000]
pxi2 = [ 0.100, 0.200]
yi2  = [ 2.000, -1.000]
pyi2 = [ 0.200, -0.100]
zi2  = [ -0.100, 0.050]
pzi2 = [ 0.900, 0.950]

xifinal2 = [1.0030362950722e00, 2.0000267496633e00]
pxifinal2 = [1.4172200742469e-01, 2.6696618573826e-01]
yifinal2 = [2.0021538002903e00, -9.9957014160228e-01]
pyifinal2 = [2.6763162365240e-01, -1.2925815827078e-01]
zifinal2 = [-1.0016187584062e-01, 5.0261547104211e-02]
pzifinal2 = [9.0534337889019e-01, 9.5346678822124e-01]

xi5  = [ -2.0, -1.0, 0.0, 1.0, 2.0 ]
pxi5 = [ -0.2, -0.1, 0.0, 0.1, 0.2 ]
yi5  = [ -1.0, -0.5, 0.0, 0.5, 1.0 ]
pyi5 = [ -0.1, -0.05, 0.0, 0.05, 0.1 ]
zi5  = [ -0.2, -0.1, 0.0, 0.1, 0.2 ]
pzi5 = [  0.8,  0.9, 1.0, 1.1, 1.2 ]

xifinal5  = [ -1.9999928778399e00, -9.9999693422668e-01, 0.0000000000000e00, 9.9999767656960e-01, 1.9999959086451e00 ]
pxifinal5 = [ -1.9999971703503e-01, -9.9999869032728e-02, 0.0000000000000e00, 9.9999888995252e-02, 1.9999979566588e-01 ]
yifinal5  = [ -9.9999272636527e-01, -4.9999683047563e-01, 0.0000000000000e00, 4.9999758493982e-01, 9.9999575931011e-01 ]
pyifinal5 = [ -9.9999725605715e-02, -4.9999872423511e-02, 0.0000000000000e00, 4.9999891831906e-02, 9.9999801293908e-02 ]
zifinal5  = [ -1.9999880462514e-01, -9.9999755289907e-02, -4.6629919381053e-11, 1.0000016810021e-01, 2.0000056466353e-01 ]
pzifinal5 = [ 7.9999970517411e-01,  8.9999971214914e-01, 9.9999971407675e-01, 1.0999997119234e00, 1.1999997066717e00 ]

xi9  = [ -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
pxi9 = [ -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
yi9 = [ -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
pyi9 = [ -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2]
zi9  = [ 0, -0.25, -0.2, -0.15, -0.1, 0.1, 0.15, 0.2, 0.25 ]
pzi9 = [  0.7, 0, 0.8, 0.85, 0.9, 1.1, 1.15, 1.2, 1.25]


xifinal9  = [ -3.9999753385049e00, -2.9999569379218e00, -1.9999647857388e00, -9.9995632170188e-01, 0.0000000000000e00, 9.9996119481153e-01, 1.9999693941762e00, 2.9999755791554e00, 3.9999799149376e00]
pxifinal9 = [ -4.0000269425828e-01, -3.0000270703508e-01, -2.0000425563770e-01, -1.0000544317117e-01, 0.0000000000000e00, 1.0000544385783e-01, 2.0000437463523e-01, 3.0000353983791e-01, 4.0000293910803e-01]
yifinal9  = [ -1.9999915140317e00, -1.4999849101804e00, -9.9999234551839e-01, -4.9999455024194e-01, 0.0000000000000e00, 4.9999536784255e-01, 9.9999375470356e-01, 1.4999934179950e00, 1.9999935925840e00]
pyifinal9 = [ -2.0000091992571e-01, -1.5000093348572e-01, -1.0000093057246e-01, -5.0000689681782e-02, 0.0000000000000e00, 5.0000659902514e-02, 1.0000089947270e-01, 1.5000095574797e-01, 2.0000093531613e-01]
zifinal9  = [ 6.8010754934456e-06, -2.4998481784087e-01, -1.9999566201089e-01, -1.4999749165298e-01, -1.0000000000001e-01, 1.0000195820814e-01, 1.5000313756545e-01, 2.0000377890961e-01, 2.5000414024036e-01]
pzifinal9 = [  7.0000053459203e-01, 6.3787411781476e-07, 8.0000064784051e-01, 8.5000076817115e-01, 9.0000095905112e-01, 1.1000007621174e00, 1.1500006333232e00, 1.2000005527379e00, 1.2500004990846e00]

@testset "BeamBeamTracking" begin
    @testset "Particles" begin
        v = [ xi pxi yi pyi zi pzi ]
        bunch = Bunch(v)
        BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (1.0e8, 1.0e11, 
                                                    1.0, 0.1, 0.1, 0.0, 1.0e14, 1, 0.0, 0.0, 0.0, 0.0, 
                                                    100.0, 0.0, 1.0, 0.0, 1.0, 0.0, [0, 0.0, 0.0, 0.0, 0.0], 0.0)))
        @test v[:,BeamTracking.XI]  ≈  xifinal
        @test v[:,BeamTracking.YI]  ≈  yifinal
        @test v[:,BeamTracking.ZI]  ≈  zifinal
        @test v[:,BeamTracking.PXI] ≈ pxifinal
        @test v[:,BeamTracking.PYI] ≈ pyifinal
        @test v[:,BeamTracking.PZI] ≈ pzifinal
        
        v = [ xi2 pxi2 yi2 pyi2 zi2 pzi2 ]
        bunch = Bunch(v)
        BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (5.0e7, 7.0e9, 
                                                    1.0, 0.1, 0.2, 2.0, 3.0e15, 2, 0.0, 0.0, 0.0, 0.0, 
                                                    0.0, 0.0, 0.4, 1.2, 0.3, 1.1, [0, 0.0, 0.0, 0.0, 0.0], 0.0)))
        @test v[:,BeamTracking.XI]  ≈  xifinal2
        @test v[:,BeamTracking.YI]  ≈  yifinal2
        @test v[:,BeamTracking.ZI]  ≈  zifinal2
        @test v[:,BeamTracking.PXI] ≈ pxifinal2
        @test v[:,BeamTracking.PYI] ≈ pyifinal2
        @test v[:,BeamTracking.PZI] ≈ pzifinal2

        v = [ xi5 pxi5 yi5 pyi5 zi5 pzi5 ]
        bunch = Bunch(v)
        BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (1.0e8, 4.0e11, 
                                                    1.0, 0.3, 0.2, 3.0, 1.0e12, 5, 0.0, 0.0, 0.0, 0.0, 
                                                    -100.0, 0.0, 1.1, 0.5, 1.3, 0.2, [0, 0.0, 0.0, 0.0, 0.0], 0.0)))
        @test v[:,BeamTracking.XI]  ≈ xifinal5
        @test v[:,BeamTracking.YI]  ≈ yifinal5
        @test v[:,BeamTracking.ZI][1:2] ≈ zifinal5[1:2]
        @test v[:,BeamTracking.ZI][3] ≈ zifinal5[3] (atol = 1e-14) #numbers too small for ≈
        @test v[:,BeamTracking.ZI][3:5] ≈ zifinal5[3:5]
        @test v[:,BeamTracking.PXI] ≈ pxifinal5
        @test v[:,BeamTracking.PYI] ≈ pyifinal5
        @test v[:,BeamTracking.PZI] ≈ pzifinal5

        v = [ xi9 pxi9 yi9 pyi9 zi9 pzi9 ]
        bunch = Bunch(v)
        BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (1.0e10, 1.0e10, 
                                                    1.0, 0.1, 0.4, 2.2, 1.0e14, 6, 0.0, 0.0, 0.0, 0.0, 
                                                    30.0, 0.0, 2.1, 0.1, 1.2, 0.0, [0, 0.0, 0.0, 0.0, 0.0], 0.0)))
        @test v[:,BeamTracking.XI] ≈ xifinal9
        @test v[:,BeamTracking.YI]  ≈ yifinal9
        @test v[:,BeamTracking.ZI] ≈ zifinal9
        @test v[:,BeamTracking.PXI] ≈ pxifinal9
        @test v[:,BeamTracking.PYI] ≈ pyifinal9
        @test v[:,BeamTracking.PZI] ≈ pzifinal9
    end
end

#Future tests with TPS and central finite differences

# xifinal = [9.95645e-1]
# pxifinal = [1.00168e-1]
# yifinal = [-4.80134e-3]
# pyifinal = [1.84375e-4]
# zifinal = [2.29639e-4]
# pzifinal = [9.00032e-1]
# @testset "BeamBeamTracking" begin
#     @testset "BeamBeam" begin
#         v = [ xi pxi yi pyi zi pzi ]
#         bunch = Bunch(v)
#         BeamTracking.launch!(bunch.coords, KernelCall(BeamTracking.track_beambeam!, (1.0e8, 1.0e11, 
#                                                     1.0, 0.1, 0.1, 0.0, 1.0e14, 1, 1.0, 2.0, 0.1, 0.2, 
#                                                     100.0, 0.5, 1.0, 0.0, 1.0, 0.0, [1, 0.1, 0.01, 0.2, 0.3], 0.6)))
#         print(bunch.coords)
#         @test v[:,BeamTracking.XI]  ≈  xifinal (rtol=5.e-3)
#         @test v[:,BeamTracking.YI]  ≈  yifinal (rtol=5.e-3)
#         @test v[:,BeamTracking.ZI]  ≈  zifinal (rtol=5.e-3)
#         @test v[:,BeamTracking.PXI] ≈ pxifinal (rtol=5.e-3)
#         @test v[:,BeamTracking.PYI] ≈ pyifinal (rtol=5.e-3)
#         @test v[:,BeamTracking.PZI] ≈ pzifinal (rtol=5.e-3)
#     end
# end
# bmad_Standard = [
#   9.8560530586828e-01  -7.1254839096554e-01   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   2.8793087232209e-04   1.0143967665022e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   9.8560530586828e-01  -7.1254839096554e-01   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   2.8793087232209e-04   1.0143967665022e+00   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   1.0000000000000e+00   5.6370935778527e-07;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   1.0000000000000e+00
# ]

# symp_Lie_PTC = [
#  -1.3399457480472e+01  -1.4399457480472e+03   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   1.4399457480472e-01   1.5399457480472e+01   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00  -1.3399457480472e+01  -1.4399457480472e+03   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   1.4399457480472e-01   1.5399457480472e+01   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   1.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   1.0000000000000e+00
# ]

# tracking = [
#   9.9712168133282e-01  -1.3672070297984e-01   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   5.7573769687519e-05   1.0028787330522e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   9.9712168133282e-01  -1.3672070297984e-01   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   5.7573769687519e-05   1.0028787330522e+00   0.0000000000000e+00   0.0000000000000e+00;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   1.0000000018792e+00   5.6400202828107e-07;
#   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   0.0000000000000e+00   5.7560564387799e-10   1.0000000000015e+00
# ]

# tps_old = [9.9712203772556085e-01 -1.3669441753244627e-01 0.0 0.0 0.0 0.0;
#        5.7570205995608794e-05 1.0028783766079528e+00 0.0 0.0 0.0 0.0;
#        0.0 0.0 9.9712203772556085e-01 -1.3669441753244627e-01 0.0 0.0;
#        0.0 0.0 5.7570205995608794e-05 1.0028783766079528e+00 0.0 0.0;
#        0.0 0.0 0.0 0.0 1.0000000000016185e+00 9.3950165115153951e-07;
#        0.0 0.0 0.0 0.0 -1.1510848234131111e-09 1.0000000000015028e+00]

# tps = [9.9712168133475609e-01 -1.3672070286721549e-01 0.0 0.0 0.0 0.0;
#        5.7577333743205789e-05 1.0028788221167801e+00 0.0 0.0 0.0 0.0;
#        0.0 0.0 9.9712159228085262e-01 -1.3672092225046839e-01 0.0 0.0;
#        0.0 0.0 5.7577333743205789e-05 1.0028788221167801e+00 0.0 0.0;
#        0.0 0.0 0.0 0.0 1.0000000075181319e+00 5.6370171554165870e-07;
#        0.0 0.0 0.0 0.0 -1.1512986031609689e-09 9.9999999999849687e-01]

# # Define a function that takes two 6x6 matrices and returns their element-wise difference
# function matrix_relative_difference(A::Matrix{<:Number}, B::Matrix{<:Number}; percent::Bool=false)
#     @assert size(A) == size(B) "Matrices must have the same size"

#     # Avoid divide-by-zero by using broadcasting and conditional replacement
#     # rel_diff = (B .- A) ./ A
#     rel_diff = abs.((B .- A))
#     rel_diff[isnan.(rel_diff)] .= 0.0        # handle 0/0 cases
#     rel_diff[isinf.(rel_diff)] .= 0.0        # handle division by zero

#     return percent ? rel_diff .* 100 : rel_diff
# end
# function pretty_print_matrix(A::AbstractMatrix; digits::Int=6, label::String="")
#     if !isempty(label)
#         println("\n$label:")
#     end

#     rows, cols = size(A)

#     # Convert each element to a scientific-notation string manually
#     formatted = Matrix{String}(undef, rows, cols)
#     for i in 1:rows, j in 1:cols
#         formatted[i, j] = lowercase(string(round(A[i, j], sigdigits=digits)))
#     end

#     # Compute the widest string in each column
#     colwidths = [maximum(length.(formatted[:, j])) for j in 1:cols]

#     # Print each row with padding
#     for i in 1:rows
#         for j in 1:cols
#             print(rpad(formatted[i, j], colwidths[j] + 2))
#         end
#         println()
#     end
# end

# function print_jacobian(J; digits=15, label="Jacobian")
#     println("\n$label:")
#     rows, cols = size(J)

#     # convert values to scientific-notation strings
#     formatted = Matrix{String}(undef, rows, cols)
#     for i in 1:rows, j in 1:cols
#         v = J[i, j]
#         if v == 0
#             formatted[i, j] = "0.0"
#         else
#             exp = floor(Int, log10(abs(v)))
#             base = round(v / 10.0^exp, sigdigits=digits)
#             formatted[i, j] = string(base, "e", exp >= 0 ? "+" : "", exp)
#         end
#     end

#     # compute column widths
#     colwidths = [maximum(length.(formatted[:, j])) for j in 1:cols]

#     # print aligned matrix
#     for i in 1:rows
#         for j in 1:cols
#             print(rpad(formatted[i, j], colwidths[j] + 2))
#         end
#         println()
#     end
# end
# function beambeam_map(x)
#     b = Bunch(copy(x)) 
#     BeamTracking.launch!(b.coords, KernelCall(
#         BeamTracking.track_beambeam_brent!,
#         (1.0e8, 1.0e11,
#          1.0, 0.1, 0.1, 1.0, 1.0e14, 1,
#          0.0, 0.0, 0.0, 0.0,
#          100.0, 0.0, 1.0, 0.0, 1.0, 0.0,
#          [0, 0.0, 0.0, 0.0, 0.0], 0.0)
#     ))
#     return b.coords 
# end

# d_z = Descriptor(6, 1)
# v5 = [0.001 -0.001 0.003 0.0002 2.0 0.6]   



# @testset "BeamBeamTracking" begin
#     @testset "BeamBeam" begin
#         v = [ xi pxi yi pyi zi pzi ]
#         bunch = Bunch(v)

#         x0 = vec(copy(v5))  # initial coordinate vector
#         fdm = central_fdm(3, 1) 
#         J = jacobian(fdm, beambeam_map, x0)


#         b1 = Bunch(v5 + transpose(@vars(d_z)))

#         BeamTracking.launch!(b1.coords, KernelCall(BeamTracking.track_beambeam!, (1.0e8, 1.0e11,
#          1.0, 0.1, 0.1, 1.0, 1.0e14, 1,
#          0.0, 0.0, 0.0, 0.0,
#          100.0, 0.0, 1.0, 0.0, 1.0, 0.0,
#          [0, 0.0, 0.0, 0.0, 0.0], 0.0)))


#         print_jacobian(J[1][2:end, :])
#         print_jacobian(GTPSA.jacobian(b1.coords.v); label="TPS Jacobian")
#         pretty_print_matrix(matrix_relative_difference(GTPSA.jacobian(b1.coords.v),J[1][2:end, :]; percent=false), digits=15, label="Absolute Difference between FiniteDiff and TPS")
#         @test v[:,BeamTracking.XI]  ≈  xifinal (rtol=5.e-4)
#         @test v[:,BeamTracking.YI]  ≈  yifinal (rtol=5.e-4)
#         @test v[:,BeamTracking.ZI]  ≈  zifinal (rtol=5.e-4)
#         @test v[:,BeamTracking.PXI] ≈ pxifinal (rtol=5.e-4)
#         @test v[:,BeamTracking.PYI] ≈ pyifinal (rtol=5.e-4)
#         @test v[:,BeamTracking.PZI] ≈ pzifinal (rtol=5.e-4)
#     end
# end

