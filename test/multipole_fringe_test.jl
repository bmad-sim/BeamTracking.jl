using LinearAlgebra: norm

@testset "Multipole hard edge fringe" begin
  J = [ 0 1  0 0  0 0;
       -1 0  0 0  0 0;
        0 0  0 1  0 0;
        0 0 -1 0  0 0;
        0 0  0 0  0 1;
        0 0  0 0 -1 0]

  # Independent, gauge free check of the generator: first order momentum kick from the
  # Lorentz force along a straight line through the edge of the Maxwell consistent field
  #   ψ = Σ_j c_j S^(2j)(s) r^(2j) Im(K w^m)/m!,  ∫ S^(k)(s) g(s) ds -> (-1)^(k-1) g^(k-1)(0)
  # must equal (-∂W/∂x, -∂W/∂y) through fifth order in the angles.
  cj(j, m) = (-1)^j*factorial(m)/(4^j*factorial(j)*factorial(j + m))
  for m in 2:5, (kn, ks) in ((1.3, 0.0), (0.0, 0.7), (0.4, -0.9))
    d = Descriptor(4, m + 6)
    x, y, xp, yp = @vars(d)
    D(f) = xp*deriv(f, 1) + yp*deriv(f, 2)
    ψj(j) = (x^2 + y^2)^j*imag((kn + im*ks)*complex(x, y)^m)/factorial(m)
    Fx = Dict{Int,Any}()
    Fy = Dict{Int,Any}()
    for j in 0:3
      c = cj(j, m)
      ψ = ψj(j)
      Fx[2j]   = get(Fx, 2j, 0) - c*deriv(ψ, 2)
      Fy[2j]   = get(Fy, 2j, 0) + c*deriv(ψ, 1)
      Fx[2j+1] = get(Fx, 2j+1, 0) + yp*c*ψ
      Fy[2j+1] = get(Fy, 2j+1, 0) - xp*c*ψ
    end
    dpx = 0
    dpy = 0
    for k in 1:6
      gx = Fx[k]
      gy = Fy[k]
      for _ in 1:k-1
        gx = D(gx)
        gy = D(gy)
      end
      dpx += (-1)^(k-1)*gx
      dpy += (-1)^(k-1)*gy
    end
    Wx, Wy, Wxp, Wyp = BeamTracking.multipole_fringe_gradient(m, kn, ks, x, y, xp, yp)
    @test normTPS(Wx + dpx) < 1e-13
    @test normTPS(Wy + dpy) < 1e-13
    # Integrability: the four components are the gradient of one function
    @test normTPS(deriv(Wx, 2) - deriv(Wy, 1)) < 1e-13
    @test normTPS(deriv(Wx, 3) - deriv(Wxp, 1)) < 1e-13
    @test normTPS(deriv(Wy, 4) - deriv(Wyp, 2)) < 1e-13
    @test normTPS(deriv(Wxp, 4) - deriv(Wyp, 3)) < 1e-13
  end

  p0 = [0.01, -0.02, 0.015, 0.01, 0.003, 0.03]
  mm = SA[2, 3, 4]
  kn = SA[0.8, 12.0, -150.0]
  ks = SA[0.3, -5.0, 90.0]

  function track_fringe(v0, mm, kn, ks, sign)
    v = reshape(copy(v0), 1, 6)
    coords = Coords(fill(STATE_ALIVE, 1), v, nothing, nothing, ())
    BeamTracking.launch!(coords, make_kernel_call(BeamTracking.multipole_fringe!, (mm, kn, ks, sign)); use_KA=false)
    return vec(coords.v)
  end

  for sign in (1, -1)
    # Symplectic
    v = transpose(p0 .+ @vars(D1))
    coords = Coords(fill(STATE_ALIVE, 1), v, nothing, nothing, ())
    BeamTracking.launch!(coords, make_kernel_call(BeamTracking.multipole_fringe!, (mm, kn, ks, sign)))
    M = GTPSA.jacobian(coords.v)[1:6,1:6]
    @test norm(transpose(M)*J*M - J) < 1e-14

    # Float64 tracking matches the TPS constant part
    @test track_fringe(p0, mm, kn, ks, sign) ≈ scalar.(vec(coords.v))

    # SIMD
    N = 37
    vs = 0.02 .* (rand(N, 6) .- 0.5)
    ref = reduce(vcat, permutedims(track_fringe(vs[j,:], mm, kn, ks, sign)) for j in 1:N)
    cs = Coords(fill(STATE_ALIVE, N), copy(vs), nothing, nothing, ())
    BeamTracking.launch!(cs, make_kernel_call(BeamTracking.multipole_fringe!, (mm, kn, ks, sign)); use_KA=false)
    @test cs.v ≈ ref

    @test_opt BeamTracking.multipole_fringe!(1, Coords(fill(STATE_ALIVE, 1), reshape(copy(p0), 1, 6), nothing, nothing, ()), mm, kn, ks, sign)
  end

  # Entrance followed by exit is the identity to first order in the strengths
  vrt = track_fringe(track_fringe(p0, mm, kn, ks, 1), mm, kn, ks, -1)
  @test maximum(abs, vrt .- p0) < 1e-10

  # Quadrupole at zero angles: the position shift equals the Lee-Whiting/Forest hard edge map
  k1 = 0.8
  for sign in (1, -1)
    v0 = [0.01, 0.0, 0.007, 0.0, 0.0, 0.03]
    rel_p = 1 + v0[6]
    x, y = v0[1], v0[3]
    vf = track_fringe(v0, SA[2], SA[k1], SA[0.0], sign)
    @test vf[1] ≈ x + sign*k1/(12*rel_p)*(x^3 + 3*y^2*x)
    @test vf[3] ≈ y - sign*k1/(12*rel_p)*(y^3 + 3*x^2*y)
  end

  # Skew quadrupole equals a rotated normal quadrupole
  rot(v, t) = [v[1]*cos(t) + v[3]*sin(t), v[2]*cos(t) + v[4]*sin(t),
              -v[1]*sin(t) + v[3]*cos(t), -v[2]*sin(t) + v[4]*cos(t), v[5], v[6]]
  t = -pi/4  # the skew quad field at w is the normal quad field at exp(i*pi/4)*w
  @test track_fringe(p0, SA[2], SA[0.0], SA[k1], 1) ≈ rot(track_fringe(rot(p0, t), SA[2], SA[k1], SA[0.0], 1), -t)

  # Differentiable with respect to the strengths, including at zero strength
  f(k) = track_fringe(eltype(k).(p0), SA[2, 3], SA[k[1], k[2]], SA[k[3], k[4]], 1)
  for k0 in ([0.0, 0.0, 0.0, 0.0], [0.8, 12.0, 0.3, -5.0])
    Jad = ForwardDiff.jacobian(f, k0)
    h = 1e-6
    Jfd = reduce(hcat, [(f(k0 .+ h .* (1:4 .== j)) .- f(k0 .- h .* (1:4 .== j))) ./ (2h) for j in 1:4])
    @test all(isfinite, Jad)
    @test maximum(abs, Jad .- Jfd) < 1e-9
  end

  # Through Beamlines: sextupole and higher fringes only with multipole_fringe_on = true
  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 5e9)
  for TM in (DriftKick, Symplectic)
    M = map(((Fringe.BothEnds, true), (Fringe.BothEnds, false), (Fringe.NoEnd, false))) do (fa, on)
      ele = LineElement(L=0.5, Kn2=20.0, Ks3=-300.0, tracking_method=TM(order=4, n_steps=4, fringe_at=fa, multipole_fringe_on=on))
      bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
      b0 = Bunch(transpose(p0 .+ @vars(D1)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
      track!(b0, bl)
      GTPSA.jacobian(b0.coords.v)[1:6,1:6]
    end
    @test norm(transpose(M[1])*J*M[1] - J) < 1e-12
    @test norm(M[1] - M[3]) > 1e-6
    @test M[2] == M[3]
  end

  # A quadrupole always gets the fringe; its sextupole component only when switched on
  M = map((true, false)) do on
    ele = LineElement(L=0.5, Kn1=0.5, Kn2=20.0, tracking_method=MatrixKick(order=4, multipole_fringe_on=on))
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=Species("electron"))
    b0 = Bunch(transpose(p0 .+ @vars(D1)), p_over_q_ref=p_over_q_ref, species=Species("electron"))
    track!(b0, bl)
    GTPSA.jacobian(b0.coords.v)[1:6,1:6]
  end
  @test norm(transpose(M[1])*J*M[1] - J) < 1e-12
  @test norm(transpose(M[2])*J*M[2] - J) < 1e-12
  @test norm(M[1] - M[2]) > 1e-6
end
