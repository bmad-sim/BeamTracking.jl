@testset "SplitIntegration" begin
  @testset "Kernels" begin
    function multipole_args(::Type{T}) where {T}
      ms = [1, 6]
      Kn0L = T(1)
      Ks0L = T(0)
      Kn5L = T(1000)
      Ks5L = T(0)
      knl = [Kn0L, Kn5L]
      ksl = [Ks0L, Ks5L]
      return ms, knl, ksl, 1, 1
    end
    
    function sk_args(::Type{T}) where {T}
      L = T(2)
      Ksol = T(0.1)
      Kn1 = T(0.1) 
      Ks1 = T(0)
      mm = [2]
      kn = [Kn1]
      sn = [Ks1]
      p0c = T(10e6)
      mc2 = T(ELECTRON.mass)
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, Ksol, mm, kn, sn, L
    end

    function mk_args(::Type{T}) where {T}
      L = T(2)
      k1 = T(0.1) 
      tilt = T(pi/4)
      w = ExactTracking.w_matrix(0,0,-tilt)
      w_inv = ExactTracking.w_inv_matrix(0,0,-tilt)
      mm = [2]
      kn = [k1*cos(-2*tilt)]
      ks = [k1*sin(-2*tilt)]
      p0c = T(10e6)
      mc2 = T(ELECTRON.mass)
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, w, w_inv, k1, mm, kn, ks, L
    end

    function dk_args(::Type{T}) where {T}
      L = T(2)
      mm = [3, 5]
      Kn2 = T(10)
      Ks2 = T(0)
      Kn4 = T(100)
      Ks4 = T(0)
      kn = [Kn2, Kn4]
      ks = [Ks2, Ks4]
      p0c = T(10e6)
      mc2 = T(ELECTRON.mass)
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, mm, kn, ks, L
    end

    function integrator_args(::Type{T}) where {T}
      L = T(2)
      mm = [3, 5]
      Kn2 = T(10)
      Ks2 = T(0)
      Kn4 = T(100)
      Ks4 = T(0)
      kn = [Kn2, Kn4]
      ks = [Ks2, Ks4]
      p0c = T(10e6)
      mc2 = T(ELECTRON.mass)
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      params = (beta_0, gamsqr_0, tilde_m, mm, kn, ks)
      ker = IntegrationTracking.dkd_multipole!
      num_steps = 1
      ds_step = T(2)
      return ker, params, ds_step, num_steps, L
    end
    
    # Scalar parameters
    test_map("bmad_maps/thin_dipole.jl", KernelCall(ExactTracking.multipole_kick!, multipole_args(Float64)); tol=1e-14)
    test_map("bmad_maps/sol_quad.jl", KernelCall(IntegrationTracking.sks_multipole!, sk_args(Float64)); tol=5e-10)
    test_map("bmad_maps/skew_quad_mk.jl", KernelCall(IntegrationTracking.mkm_quadrupole!, mk_args(Float64)); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.dkd_multipole!, dk_args(Float64)); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.order_two_integrator!, integrator_args(Float64)); tol=5e-10)
    test_map("bmad_maps/order_four.jl", KernelCall(IntegrationTracking.order_four_integrator!, integrator_args(Float64)); tol=5e-10)
    test_map("bmad_maps/order_six.jl", KernelCall(IntegrationTracking.order_six_integrator!, integrator_args(Float64)); tol=5e-9)
    test_map("bmad_maps/order_eight.jl", KernelCall(IntegrationTracking.order_eight_integrator!, integrator_args(Float64)); tol=5e-9)

    # GTPSA parameters
    test_map("bmad_maps/thin_dipole.jl", KernelCall(ExactTracking.multipole_kick!, multipole_args(TPS64{D10})); tol=1e-14)
    test_map("bmad_maps/sol_quad.jl", KernelCall(IntegrationTracking.sks_multipole!, sk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/skew_quad_mk.jl", KernelCall(IntegrationTracking.mkm_quadrupole!, mk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.dkd_multipole!, dk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.order_two_integrator!, integrator_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/order_four.jl", KernelCall(IntegrationTracking.order_four_integrator!, integrator_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/order_six.jl", KernelCall(IntegrationTracking.order_six_integrator!, integrator_args(TPS64{D10})); tol=5e-9)
    test_map("bmad_maps/order_eight.jl", KernelCall(IntegrationTracking.order_eight_integrator!, integrator_args(TPS64{D10})); tol=5e-9)
  end
end
