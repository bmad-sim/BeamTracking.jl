#using Plots
include("../ext/BeamTrackingBeamlinesExt/utils_bl.jl")


@testset "ImplicitIntegration" begin
  dist(v1, v2) = sqrt(sum(abs2, v1 - v2))

  function A_sol(k)
    func = (x, y, s, t) -> (0, -k*y/2, k*x/2, 0)
    jac = (x, y, s, t) -> SA[
      0 0 0 0;
      0 -k/2 0 0;
      k/2 0 0 0;
      0 0 0 0
    ]
    return func, jac
  end

  function A_dipole(g) 
    func = (x, y, s, t) -> (0, 0, 0, -g*x*(1+g*x/2))
    jac = (x, y, s, t) -> SA[
      0 0 0;
      0 0 0;
      0 0 0;
      -g*(1+g*x) 0 0;
    ]
    return func, jac
  end

  function A_quad(k) 
    func = (x, y, s, t) -> (0, 0, 0, -k*(x*x-y*y)/2)
    jac = (x, y, s, t) -> SA[
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      -k*x k*y 0 0
    ]
    return func, jac
  end

  @inline function k_to_Az(ms, knl, ksl, x, y)
    z = complex(x, y)
    az = zero(x * y)
    @inbounds for j in eachindex(ms, knl, ksl)
      m = ms[j]
      coeff = complex(knl[j], ksl[j]) / factorial(m)
      az -= real(coeff * z^m)
    end
    return az
  end

  @inline function k_to_∂Az(ms, knl, ksl, x, y)
    z = complex(x, y)
    sum_d = zero(complex(x, y))
    @inbounds for j in eachindex(ms, knl, ksl)
      m = ms[j]
      if m > 0
        coeff = complex(knl[j], ksl[j]) / factorial(m - 1)
        sum_d += coeff * z^(m - 1)
      end
    end
    dAz_dx = -real(sum_d)
    dAz_dy = imag(sum_d)
    return dAz_dx, dAz_dy
  end

  function _A(ms, knl, ksl)
    @assert length(ms) == length(knl) == length(ksl) "ms, knl, and ksl must have the same length"
    A = (x, y, s, t) -> (0, 0, 0, k_to_Az(ms, knl, ksl, x, y))
    return A
  end

  function _JA(ms, knl, ksl)
    @assert length(ms) == length(knl) == length(ksl) "ms, knl, and ksl must have the same length"
    JA = (x, y, s, t) -> begin
      dAz_dx, dAz_dy = k_to_∂Az(ms, knl, ksl, x, y)
      SA[
        0      0      0 0;
        0      0      0 0;
        0      0      0 0;
        dAz_dx dAz_dy 0 0
      ]
    end
    return JA
  end


  q0 = SA[0.0031, -0.0012, 2.312e-8];
  p0 = SA[0.001, -0.002, 0.001355];
  ds = 0.013131;
  k = 1.23;
  g = 0.2512;
  mm = SA[1,2,4]; kn = SA[0.25, -0.41, 0.91]; ks = SA[0.0, 0.0 ,0.0];

  species = BeamTracking.Species("proton")
  a = gyromagnetic_anomaly(species)
  p0c = 1e11
  p_over_q_ref = p0c / BeamTracking.C_LIGHT
  tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(species, p_over_q_ref)
  tilde_m2 = tilde_m^2





  order_rows = (
    BeamTracking.order_two_integrator!,
    BeamTracking.order_four_integrator!,
    BeamTracking.order_six_integrator!,
    BeamTracking.order_eight_integrator!,
  )



  sol_exact_kc(k, L) = BeamTracking.KernelCall(BeamTracking.exact_solenoid!, (k, beta_0, gamsqr_0, tilde_m, L))
  dip_exact_kc(g, L) = BeamTracking.KernelCall(BeamTracking.exact_bend!, (g * L, g, g, tilde_m, beta_0, L))



  function order_kcall(order, ker, params, L, num_steps=1)
    @assert order in [2,4,6,8]    
    BeamTracking.KernelCall(order_rows[Int(order/2)], (ker, params, nothing, L/num_steps, num_steps, nothing, Val(false), Val(false), L))
  end


  function multipole_kcall(mm,kn,ks,L; ker=BeamTracking.dkd_multipole!, order=6, num_steps=1_000)
    params = (species.charge, species.mass, false, beta_0, gamsqr_0, tilde_m, gyromagnetic_anomaly(species), mm, kn, ks)
    return order_kcall(order, ker, params, L, num_steps)
  end

  function imp_solenoid_kcall(ksol, order, impl_kernel, L, num_steps)
    A, JA = A_sol(ksol)
    args = (A, JA, [0.0], a, 0.0, beta_0, tilde_m2,
            1e-14, 1e-14, 50)
    return order_kcall(order, impl_kernel, args, L, num_steps)
  end

  function imp_dipole_kcall(g, order, impl_kernel, L, num_steps)
    A, JA = A_dipole(g)
    args = (A, JA, [0.0], a, g, beta_0, tilde_m2, 
            1e-14, 1e-14, 15)
    return order_kcall(order, impl_kernel, args, L, num_steps)
  end

  function imp_multipole_kcall(mm, kn, ks, order, impl_kernel, L, num_steps)
    A = _A(mm,kn,ks)
    JA = _JA(mm,kn,ks)
    args = (A, JA, [0.0], a, 0.0, beta_0, tilde_m2,
            1e-14, 1e-14, 15)
    return order_kcall(order, impl_kernel, args, L, num_steps)
  end

  function imp_wiggler_kcall(params, order, impl_kernel, L, num_steps)
    _A = params[1]
    k = params[2]
    r = params[3]
  
    function A(x,y,s,t) 
      shx = sinh(k[1]*(x+r[1]))
      shy = sinh(k[2]*(y+r[2])); chy = cosh(k[2]*(y+r[2]))
      sz = sin(k[3]*s+r[3]); cz = cos(k[3]*s+r[3])
      return -_A/k[1] * [
                      0,
                      0,
                      shx*shy*sz,
                      k[2]/k[3]*shx*chy*cz
      ]
    end

    function JA(x,y,s,t) 
      shx = sinh(k[1]*(x+r[1])); chx = cosh(k[1]*(x+r[1]))
      shy = sinh(k[2]*(y+r[2])); chy = cosh(k[2]*(y+r[2]))
      sz = sin(k[3]*s+r[3]); cz = cos(k[3]*s+r[3])
      return -_A/k[1] * [
        0 0 0 0
        0 0 0 0
        k[1]*chx*shy*sz k[2]*shx*chy*sz 0 0
        k[1]*k[2]/k[3]*chx*chy*cz k[2]^2/k[3]*shx*shy*cz 0 0
      ]
    end

    args = (A, JA, [r[4]], 0.0, beta_0, tilde_m2,
        1e-14, 1e-14, 15)
    return order_kcall(order, impl_kernel, args, L, num_steps)
  end

  magnets = Dict(
    :Solenoid => ((L)->sol_exact_kc(k, L), (args...)->imp_solenoid_kcall(k, args...)),
    :Dipole => ((L)->dip_exact_kc(g, L), (args...)->imp_dipole_kcall(g, args...)),
    :Multipole => ((L)->multipole_kcall(mm,kn,ks,L), (args...)->imp_multipole_kcall(mm, kn, ks, args...)),
    :Wiggler => ((L)->error("No reference defined for wiggler"), (args...) -> imp_wiggler_kcall([
        8789.313/b.p_over_q_ref,
        [26.990786,0.00737,sqrt(26.990786^2+0.00737^2)],
        [0.0, 0.0, 1.570796, 0.0]
        ], args...)
    )
  )
  
  
  # =
  @testset "Vector Map" begin

    function make_bunch(q, p)
      v = reshape([q[1], p[1], q[2], p[2], q[3], p[3]], 1, 6)
      quat = reshape([1.0, 0.0, 0.0, 0.0], 1, 4)
      return Bunch(deepcopy(v), deepcopy(quat), species=species, p_over_q_ref=p_over_q_ref)
    end

    for magnet in [:Solenoid, :Dipole, :Multipole] #keys(magnets)
      @testset "$(magnet)" begin
        kc_ref_fn, kc_impl_fn = magnets[magnet]

        b_ref = make_bunch(q0, p0)
        BeamTracking.launch!(b_ref.coords, kc_ref_fn(1.0))

        for order in [2,4,6,8]
          @testset "Order $(order)" begin
    
            # plt = plot()
            # plt_q = plot()

            # _tmp = []
            # _tmp_q = []
            
            # for log_num_steps in 1:10

              num_steps = 2 .^ (18-order) #log_num_steps

              b_impl = make_bunch(q0, p0)
              kc_impl = kc_impl_fn(order, BeamTracking.symplectic_step!, 1.0, num_steps)
              BeamTracking.launch!(b_impl.coords, kc_impl)  
              # push!(_tmp, sqrt(sum((b_ref.coords.v .- b_impl.coords.v).^2)/6))
              # push!(_tmp_q, sqrt(sum((b_ref.coords.q .- b_impl.coords.q).^2)/4))

              # Orbit check
              @test dist(b_ref.coords.v, b_impl.coords.v) ≈ 0 atol=1e-11
              
              # Spin check
              if magnet == :Multipole
                @test dist(b_ref.coords.q, b_impl.coords.q) ≈ 0 atol=2.5e-10
              end

            # end

            # plot!(plt, 2 .^ collect(1:10), _tmp; 
            #       label="test", lw=2, marker=:circle, markersize=4,
            #       xscale=:log10, yscale=:log10, xlabel="num_steps", ylabel="||x - x_ref||₂", 
            #       title="$(magnet): implicit-$(order) vs reference"
            #       )

            # plot!(plt_q, 2 .^ collect(1:10), _tmp_q; 
            #       label="test", lw=2, marker=:circle, markersize=4,
            #       xscale=:log10, yscale=:log10, xlabel="num_steps", ylabel="||q - q_ref||₂", 
            #       title="$(magnet) spin: implicit-$(order) vs reference"
            #       )

            # savefig(plt, "plots_2/implicit_$(magnet)-$(order).png")
            # savefig(plt_q, "plots_2/spin_implicit_$(magnet)-$(order).png")

          end
        end
      end
    end
  end

  # = #
  @testset "TPSA Map" begin
    D_map = Descriptor(6, 3)

    function make_tps_coords(q, p, D)
      dz = @vars(D)
      v = transpose([q[1] + dz[1], p[1] + dz[2],
                     q[2] + dz[3], p[2] + dz[4],
                     q[3] + dz[5], p[3] + dz[6]])
      state = UInt8[BeamTracking.STATE_ALIVE]
      return BeamTracking.Coords(state, v, nothing, nothing)
    end
    
  
    for magnet in [:Solenoid, :Dipole, :Multipole] #keys(magnets)
      @testset "$(magnet)" begin

        kc_ref_fn, kc_impl_fn = magnets[magnet]
        c_ref = make_tps_coords(q0, p0, D_map)
        BeamTracking.launch!(c_ref, kc_ref_fn(1.0))

          for order in [2,4,6,8]
            @testset "Order $(order)" begin
              #plt = plot()
              #_tmp = []
            
              #for log_num_steps in 1:10
                num_steps = 2 .^ (16-order) #log_num_steps

                c_imp = make_tps_coords(q0, p0, D_map)
                kc_imp = kc_impl_fn(order, BeamTracking.symplectic_step_tpsa!, 1.0, num_steps)
                BeamTracking.launch!(c_imp, kc_imp)


                # n = GTPSA.numcoefs(c_ref.v[1,1])
                # err = 0
                # for i in 1:length(c_ref.v[1,:])
                #     for j in 0:n-1
                #         c1, c2 = c_ref.v[1,i][j], c_imp.v[1,i][j]
                #         err += (c1 - c2)^2
                #     end
                # end
                # err = sqrt(err/n)
                #err = sum((GTPSA.scalar.(c_ref.v[1,:]) .- GTPSA.scalar.(c_imp.v[1,:])).^2)
                #err += norm(GTPSA.jacobian(c_ref.v[1,:]) .- GTPSA.jacobian(c_imp.v[1,:]))^2
                #err += sum(norm.([GTPSA.hessian(c_ref.v[1,i]) .- GTPSA.hessian(c_imp.v[1,i]) for i in 1:6]).^2)
                #push!(_tmp, sqrt(err))

                @test coeffs_approx_equal(c_ref.v[1,:], c_imp.v[1,:], 2.5e-9)
              #end

            # plot!(plt, 2 .^ collect(1:10), _tmp; 
            #       label="test", lw=2, marker=:circle, markersize=4,
            #       xscale=:log10, yscale=:log10, xlabel="num_steps", ylabel="||x - x_ref||₂", 
            #       title="$(magnet): implicitTPS-$(order) vs reference"
            #       )

            # savefig(plt, "plots/implicitTPS-$(order)_$(magnet).png")
            end
          end
      end
    end
  end #
end # =#

