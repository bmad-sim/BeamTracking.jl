@testset "SplitIntegration" begin
  @testset "Kernels" begin
    function multipole_args(::Type{T}) where {T}
      ms = SA[1, 6]
      Kn0L = T(0.1)
      Ks0L = T(0)
      Kn5L = T(1000)
      Ks5L = T(0)
      knl = SA[Kn0L, Kn5L]
      ksl = SA[Ks0L, Ks5L]
      return ms, knl, ksl, -1
    end
    
    function sk_args(::Type{T}) where {T}
      L = T(2)
      Ksol = T(0.1)
      Kn1 = T(0.1) 
      Ks1 = T(0)
      mm = SA[2]
      kn = SA[Kn1]
      sn = SA[Ks1]
      a = T(0.00115965218046)
      p0c = T(10e6)
      mc2 = T(BeamTracking.massof(Species("electron")))
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, sn, L
    end

    function mk_args(::Type{T}) where {T}
      L = T(2)
      k1 = T(0.1) 
      tilt = T(pi/4)
      w = ExactTracking.w_quaternion(0,0,-tilt)
      w_inv = ExactTracking.w_inv_quaternion(0,0,-tilt)
      mm = SA[2]
      kn = SA[k1*cos(-2*tilt)]
      ks = SA[k1*sin(-2*tilt)]
      a = T(0.00115965218046)
      p0c = T(10e6)
      mc2 = T(BeamTracking.massof(Species("electron")))
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L
    end

    function dk_args(::Type{T}) where {T}
      L = T(2)
      mm = SA[3, 5]
      Kn2 = T(10)
      Ks2 = T(0)
      Kn4 = T(100)
      Ks4 = T(0)
      kn = SA[Kn2, Kn4]
      ks = SA[Ks2, Ks4]
      a = T(0.00115965218046)
      p0c = T(10e6)
      mc2 = T(BeamTracking.massof(Species("electron")))
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      return beta_0, gamsqr_0, tilde_m, a, mm, kn, ks, L
    end

    function bk_straight_args(::Type{T}) where {T}
      L = T(2)
      mm = SA[1, 2]
      Kn0 = T(0.1)
      Ks0 = T(0)
      Kn1 = T(0.1)
      Ks1 = T(0)
      kn = SA[Kn0, Kn1]
      ks = SA[Ks0, Ks1]
      w = w_inv = SA[1.0 0.0 0.0 0.0]
      a = T(0.00115965218046)
      p0c = T(10e6)
      mc2 = T(BeamTracking.massof(Species("electron")))
      tilde_m = mc2/p0c
      beta_0 = 1/sqrt(1 + tilde_m^2)
      params = (tilde_m, beta_0, a, 0, 0, 0, w, w_inv, Kn0, mm, kn, ks)
      ker = IntegrationTracking.bkb_multipole!
      num_steps = 10
      ds_step = T(0.2)
      return ker, params, ds_step, num_steps, L
    end

    function integrator_args(::Type{T}) where {T}
      L = T(2)
      mm = SA[3, 5]
      Kn2 = T(10)
      Ks2 = T(0)
      Kn4 = T(100)
      Ks4 = T(0)
      kn = SA[Kn2, Kn4]
      ks = SA[Ks2, Ks4]
      a = T(0.00115965218046)
      p0c = T(10e6)
      mc2 = T(BeamTracking.massof(Species("electron")))
      tilde_m = mc2/p0c
      gamsqr_0 = 1 + 1/tilde_m^2
      beta_0 = 1/sqrt(1 + tilde_m^2)
      params = (beta_0, gamsqr_0, tilde_m, a, mm, kn, ks)
      ker = IntegrationTracking.dkd_multipole!
      num_steps = 1
      ds_step = T(2)
      return ker, params, ds_step, num_steps, L
    end
    
    # Scalar parameters
    test_map("bmad_maps/thin_dipole.jl", KernelCall(ExactTracking.multipole_kick!, multipole_args(Float64)); tol=1e-14, no_scalar_allocs=true)
    test_map("bmad_maps/sol_quad.jl", KernelCall(IntegrationTracking.sks_multipole!, sk_args(Float64)); tol=5e-10, no_scalar_allocs=true)
    test_map("bmad_maps/skew_quad_mk.jl", KernelCall(IntegrationTracking.mkm_quadrupole!, mk_args(Float64)); tol=5e-10, no_scalar_allocs=true)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.dkd_multipole!, dk_args(Float64)); tol=5e-10, no_scalar_allocs=true)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.order_two_integrator!, integrator_args(Float64)); tol=5e-10, no_scalar_allocs=true)
    test_map("bmad_maps/order_four.jl", KernelCall(IntegrationTracking.order_four_integrator!, integrator_args(Float64)); tol=5e-10, no_scalar_allocs=true)
    test_map("bmad_maps/order_six.jl", KernelCall(IntegrationTracking.order_six_integrator!, integrator_args(Float64)); tol=5e-9, no_scalar_allocs=true)
    test_map("bmad_maps/order_eight.jl", KernelCall(IntegrationTracking.order_eight_integrator!, integrator_args(Float64)); tol=5e-9, no_scalar_allocs=true)
    test_map("bmad_maps/straight_dipole_bk.jl", KernelCall(IntegrationTracking.order_six_integrator!, bk_straight_args(Float64)); tol=2e-6, no_scalar_allocs=true)

    # GTPSA parameters
    test_map("bmad_maps/thin_dipole.jl", KernelCall(ExactTracking.multipole_kick!, multipole_args(TPS64{D10})); tol=1e-14)
    test_map("bmad_maps/sol_quad.jl", KernelCall(IntegrationTracking.sks_multipole!, sk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/skew_quad_mk.jl", KernelCall(IntegrationTracking.mkm_quadrupole!, mk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.dkd_multipole!, dk_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/sex_dec.jl", KernelCall(IntegrationTracking.order_two_integrator!, integrator_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/order_four.jl", KernelCall(IntegrationTracking.order_four_integrator!, integrator_args(TPS64{D10})); tol=5e-10)
    test_map("bmad_maps/order_six.jl", KernelCall(IntegrationTracking.order_six_integrator!, integrator_args(TPS64{D10})); tol=5e-9)
    test_map("bmad_maps/order_eight.jl", KernelCall(IntegrationTracking.order_eight_integrator!, integrator_args(TPS64{D10})); tol=5e-9)
    test_map("bmad_maps/straight_dipole_bk.jl", KernelCall(IntegrationTracking.order_six_integrator!, bk_straight_args(TPS64{D10})); tol=2e-6)
  end
end
