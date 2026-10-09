using LinearAlgebra: norm

@testset "Bend second order fringe" begin
  # Hwang-Lee generator, Bmad manual "Bend Second Order Fringe Map". s = +1 entrance, -1 exit.
  function hwang_omega(v, g, k1, e, I, s)
    x, px, y, py, z, pz = v
    t = tan(e); sc = sec(e); sn = sin(e)
    d = 1 + pz
    K = g*t*(x^2 - y^2)/2 + (g^2*sc^3*(1 + sn^2)*I*y^2 + x^3*(4*k1*t - g^2*t^3)/12 + x*y^2*(-4*k1*t + g^2*t*sc^2)/4)/d
    B = g*(t^2*(x^2*px - 2*x*y*py) - sc^2*y^2*px)/(2*d)
    return K + s*B
  end

  pb(f, h) = sum(deriv(f, 2j-1)*deriv(h, 2j) - deriv(f, 2j)*deriv(h, 2j-1) for j in 1:3)

  function lie_exp(v, Om; nmax=12)
    return map(v) do zi
      out = zi
      term = zi
      for n in 1:nmax
        term = pb(Om, term)/n
        out += term
      end
      out
    end
  end

  J = [ 0 1  0 0  0 0;
       -1 0  0 0  0 0;
        0 0  0 1  0 0;
        0 0 -1 0  0 0;
        0 0  0 0  0 1;
        0 0  0 0 -1 0]

  g, k1, e1, e2, I1, I2 = 0.37, 0.8, 0.3, -0.2, 0.04, 0.03
  args = (nothing, nothing, g, k1, nothing, nothing, e1, e2, I1, I2)
  D = Descriptor(6, 5)

  for (s, e, I) in ((1, e1, I1), (-1, e2, I2))
    # Agrees with exp(:Ω:) through second order
    v0 = @vars(D)
    v = transpose(copy(v0))
    coords = Coords(fill(STATE_ALIVE, 1), v, nothing, nothing, ())
    BeamTracking.launch!(coords, make_kernel_call(BeamTracking.fringe!, (args..., s)))
    M_lie = lie_exp(v0, hwang_omega(v0, g, k1, e, I, s))
    for j in 1:6, o in 0:2
      @test normTPS(getord(coords.v[1,j] - M_lie[j], o)) < 1e-14
    end

    # Exactly symplectic away from the origin
    p0 = [0.01, -0.02, 0.015, 0.01, 0.003, 0.03]
    v = transpose(p0 .+ @vars(D1))
    coords = Coords(fill(STATE_ALIVE, 1), v, nothing, nothing, ())
    BeamTracking.launch!(coords, make_kernel_call(BeamTracking.fringe!, (args..., s)))
    M = GTPSA.jacobian(coords.v)[1:6,1:6]
    @test norm(transpose(M)*J*M - J) < 1e-14

    # Float64 tracking matches the TPS constant part
    vf = reshape(copy(p0), 1, 6)
    coords_f = Coords(fill(STATE_ALIVE, 1), vf, nothing, nothing, ())
    BeamTracking.launch!(coords_f, make_kernel_call(BeamTracking.fringe!, (args..., s)); use_KA=false)
    @test vf ≈ scalar.(coords.v)

    @test_opt BeamTracking.fringe!(1, coords_f, args..., s)
  end

  # Through Beamlines: a symplectic bend with fringes at both ends
  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 5e9)
  ele = LineElement(L=1.0, g=0.1, Kn0=0.1, Kn1=0.2, e1=0.1, e2=-0.05, edge1_int=0.02, edge2_int=0.01,
                    tracking_method=BendKick(order=4, n_steps=4))
  bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
  b0 = Bunch(transpose([0.01, -0.02, 0.015, 0.01, 0.003, 0.03] .+ @vars(D1)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
  track!(b0, bl)
  M = GTPSA.jacobian(b0.coords.v)[1:6,1:6]
  @test norm(transpose(M)*J*M - J) < 1e-12
end
