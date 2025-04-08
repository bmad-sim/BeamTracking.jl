
@testset "LinearTracking" begin
  @testset "Utility functions" begin
    # Quadrupole
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

  end

  @testset "Kernels" begin
    # define default values 
    d = Descriptor(6+10, 1) # allow 10 parameters
    k = @vars(d)[7:end] # parameters

    Ls = rand(Float64)
    Lt = rand(Float64) + k[1]

    gamma_0s = 1e5*rand(Float64)
    gamma_0t = 1e5*rand(Float64) + k[2]

    K1s = rand(Float64)
    K1t = rand(Float64) + k[3]

    # Drift =======================================================================
    M_drift(L, gamma_0) = [1.0  L    0.0  0.0  0.0  0.0;
                           0.0  1.0  0.0  0.0  0.0  0.0;
                           0.0  0.0  1.0  L    0.0  0.0;
                           0.0  0.0  0.0  1.0  0.0  0.0;
                           0.0  0.0  0.0  0.0  1.0  L/gamma_0^2;
                           0.0  0.0  0.0  0.0  0.0  1.0] 
    # Scalar parameters
    test_matrix(LinearTracking.linear_drift!, M_drift(Ls,gamma_0s), Ls, Ls/gamma_0s^2)

    # GTPSA parameters
    test_matrix(LinearTracking.linear_drift!, M_drift(Lt,gamma_0t), Lt, Lt/gamma_0t^2)

    # Quadrupole ==================================================================
    function m_quad(K1,L,gamma_0)
      M = zeros(promote_type(map(t->typeof(t), (K1,L,gamma_0))...), 6, 6)
      mx, my = LinearTracking.linear_quad_matrices(K1, L)
      M[1:2,1:2] .= mx
      M[3:4,3:4] .= my
      M[5:6,5:6] .= [1.0 L/gamma_0^2; 0.0 1.0]
      return M
    end

    # Scalar parameters
    test_matrix(LinearTracking.linear_coast_uncoupled!, m_quad(K1s, Ls, gamma_0s), LinearTracking.linear_quad_matrices(K1s, Ls)..., Ls/gamma_0s^2)

    # GTPSA parameters
    test_matrix(LinearTracking.linear_coast_uncoupled!, m_quad(K1t, Lt, gamma_0t), LinearTracking.linear_quad_matrices(K1t, Lt)..., Lt/gamma_0t^2)

  end
end
