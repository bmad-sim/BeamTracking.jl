
@testset "LinearTracking" begin
  @testset "Utility functions" begin
    # Thick Quadrupole
    mf(K1,L) = [cos(sqrt(K1)*L)            sincu(sqrt(K1)*L)*L;  
                -sqrt(K1)*sin(sqrt(K1)*L)  cos(sqrt(K1)*L)     ]
    md(K1,L) = [cosh(sqrt(K1)*L)           sinhcu(sqrt(K1)*L)*L; 
                sqrt(K1)*sinh(sqrt(K1)*L)  cosh(sqrt(K1)*L)     ]

    L = 1.2
    K1 = 0.36
  
    # Focusing
    @test all(LinearTracking.linear_quad_matrices(K1, L) .== (mf(K1,L), md(K1,L)))

    # Defocusing
    @test all(LinearTracking.linear_quad_matrices(-K1, L) .== (md(K1,L), mf(K1,L)))


    # Thin Quadrupole
    mft(K1L) = [one(K1L) zero(K1L);
                -K1L     one(K1L)  ]
    mdt(K1L) = [one(K1L) zero(K1L);
                K1L      one(K1L)  ]
    K1L = K1*L

    @test all(LinearTracking.linear_thin_quad_matrices(K1L) .== (mft(K1L), mdt(K1L)))
    @test all(LinearTracking.linear_thin_quad_matrices(-K1L) .== (mdt(K1L), mft(K1L)))
    

    # Solenoid
    M_sol_ESR = [ 0.9620267534604E+00  0.2680023374209E+01  0.1911315753266E+00 0.5324561791884E+00;
                 -0.1363095540075E-01  0.9620267534604E+00 -0.2708142959206E-02 0.1911315753266E+00;
                 -0.1911315753266E+00 -0.5324561791884E+00  0.9620267534604E+00 0.2680023374209E+01;
                  0.2708142959206E-02 -0.1911315753266E+00 -0.1363095540075E-01 0.9620267534604E+00]
    @test LinearTracking.linear_solenoid_matrix(0.142634259959, 2.75) ≈ M_sol_ESR

    # SBend
    gamma_0 = 3.4924264755852841E+04
    L = 3.8
    K0 = -1.4085135130897E-3

    M_sbend_ESR = [  0.9999856762017098E+00  0.3799981856504840E+01  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00 -0.1016944328690584E-01
                    -0.7538823207846669E-05  0.9999856762017098E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00 -0.5352325794382752E-02
                    0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  0.3800000000000000E+01  0.0000000000000000E+00  0.0000000000000000E+00
                    0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  0.0000000000000000E+00  0.0000000000000000E+00
                    0.5352325794382752E-02  0.1016944328690584E-01  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01 -0.1814037965058192E-04
                    0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01]

    M = zeros(6,6)
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, gamma_0)
    M[1:2, 1:2] = mx
    M[3:4, 3:4] = my
    M[5, 5] = 1.0
    M[6, 6] = 1.0
    M[5, 6] = r56
    M[5, 1:4] = t
    M[1:4, 6] = d
    @test M ≈ M_sbend_ESR


    L=3.8
    K0=-1.4085135130897E-1
    e1 = -2.6e1
    e2 = 2e-1

    M_crazy_bend = [  0.1461364082101138E+01  0.3621145980031785E+01  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00 -0.9928997816884232E+00
                      0.2924452887724791E-01  0.7567578276535424E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00 -0.4816940474576835E+00
                      0.0000000000000000E+00  0.0000000000000000E+00  0.3690896823132183E+00  0.3800000000000000E+01  0.0000000000000000E+00  0.0000000000000000E+00
                      0.0000000000000000E+00  0.0000000000000000E+00 -0.1554907888474444E+00  0.1108497533216086E+01  0.0000000000000000E+00  0.0000000000000000E+00
                      0.6748934931787792E+00  0.9928997816884232E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01 -0.1788540168527052E+00
                      0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01]
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, gamma_0, e1, e2)
    M[1:2, 1:2] = mx
    M[3:4, 3:4] = my
    M[5, 5] = 1.0
    M[6, 6] = 1.0
    M[5, 6] = r56
    M[5, 1:4] = t
    M[1:4, 6] = d

    @test M ≈ M_crazy_bend
  end

  @testset "Kernels" begin
    Ls = rand(Float64)
    Lt = TPS{D1}(rand(Float64))

    gamma_0s = rand(Float64)
    gamma_0t = TPS{D1}(rand(Float64))

    # Drift =======================================================================
    M_drift(L, gamma_0) = SA[1.0  L    0.0  0.0  0.0  0.0;
                             0.0  1.0  0.0  0.0  0.0  0.0;
                             0.0  0.0  1.0  L    0.0  0.0;
                             0.0  0.0  0.0  1.0  0.0  0.0;
                             0.0  0.0  0.0  0.0  1.0  L/gamma_0^2;
                             0.0  0.0  0.0  0.0  0.0  1.0] 
    # Scalar parameters
    test_matrix(M_drift(Ls,gamma_0s), KernelCall(LinearTracking.linear_drift!, (Ls, Ls/gamma_0s^2)))

    # GTPSA parameters
    test_matrix(M_drift(Lt,gamma_0t), KernelCall(LinearTracking.linear_drift!, (Lt, Lt/gamma_0t^2)))

    # Coast uncoupled  =============================================================
    function coast_uncoupled(::Type{T}) where {T}
      mx = @SArray T[1 2; 3 4]
      my = @SArray T[5 6; 7 8]
      r56 = T(9)
      d = @SArray T[10, 11, 12, 13]
      t = @SArray T[14, 15, 16, 17]
      return mx, my, r56, d, t
    end

    # Scalar parameters
    mxs, mys, r56s, ds, ts = coast_uncoupled(Float64)
    Ms = zeros(6,6)
    Ms[1:2, 1:2] = mxs
    Ms[3:4, 3:4] = mys
    Ms[5, 5] = 1.0
    Ms[6, 6] = 1.0
    Ms[5, 6] = r56s
    test_matrix(Ms, KernelCall(LinearTracking.linear_coast_uncoupled!, (mxs, mys, r56s, nothing, nothing)))
    Ms[5, 1:4] = ts
    Ms[1:4, 6] = ds
    test_matrix(Ms, KernelCall(LinearTracking.linear_coast_uncoupled!, (mxs, mys, r56s, ds, ts)))

    # GTPSA parameters
    mxt, myt, r56t, dt, tt = coast_uncoupled(TPS64{D1})
    Mt = zeros(TPS64{D1},6,6)
    Mt[1:2, 1:2] = mxt
    Mt[3:4, 3:4] = myt
    Mt[5, 5] = 1.0
    Mt[6, 6] = 1.0
    Mt[5, 6] = r56t
    test_matrix(scalar.(Mt), KernelCall(LinearTracking.linear_coast_uncoupled!, (mxt, myt, r56t, nothing, nothing)))
    Mt[5, 1:4] = tt
    Mt[1:4, 6] = dt
    test_matrix(scalar.(Mt), KernelCall(LinearTracking.linear_coast_uncoupled!, (mxt, myt, r56t, dt, tt)))

    # Coast =============================================================
    function coast(::Type{T}) where {T}
      mxy = @SArray T[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
      r56 = T(17)
      d = @SArray T[18, 19, 20, 21]
      t = @SArray T[22, 23, 24, 25]
      return mxy, r56, d, t
    end

    # Scalar parameters
    mxys, r56s, ds, ts = coast(Float64)
    Ms = zeros(6,6)
    Ms[1:4, 1:4] = mxys
    Ms[5, 5] = 1.0
    Ms[6, 6] = 1.0
    Ms[5, 6] = r56s
    test_matrix(Ms, KernelCall(LinearTracking.linear_coast!, (mxys, r56s, nothing, nothing)))
    Ms[5, 1:4] = ts
    Ms[1:4, 6] = ds
    test_matrix(Ms, KernelCall(LinearTracking.linear_coast!, (mxys, r56s, ds, ts)))

    # GTPSA parameters
    mxyt, r56t, dt, tt = coast(TPS64{D1})
    Mt = zeros(TPS64{D1},6,6)
    Mt[1:4, 1:4] = mxyt
    Mt[5, 5] = 1.0
    Mt[6, 6] = 1.0
    Mt[5, 6] = r56t
    test_matrix(scalar.(Mt), KernelCall(LinearTracking.linear_coast!, (mxyt, r56t, nothing, nothing)))
    Mt[5, 1:4] = tt
    Mt[1:4, 6] = dt
    test_matrix(scalar.(Mt), KernelCall(LinearTracking.linear_coast!, (mxyt, r56t, dt, tt)))


    # 6D =============================================================
    function sixD(::Type{T}) where {T}
      SMatrix{6,6}(collect(T, reshape(1:36, 6,6)))
    end

    # Scalar parameters
    Ms = sixD(Float64)
    test_matrix(Ms, KernelCall(LinearTracking.linear_6D!, (Ms,)))

    # GTPSA parameters
    Mt = sixD(TPS64{D1})
    test_matrix(Ms, KernelCall(LinearTracking.linear_6D!, (Mt,)))
  end
end