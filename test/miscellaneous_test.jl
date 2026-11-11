using Test
using AtomicAndPhysicalConstants
import BeamTracking as BT

@testset "Miscellaneous" begin
  @test Beamlines.R_to_v(Species("He--"), -0.3) ≈ 1.4442639106654184e7

  # Energy conversions must be consistent for any charge state (|Z| != 1 included).
  @testset "Energy conversions $name" for (name, R) in (("electron", -1.0), ("proton", 1.0),
                                    ("He++", 1.0), ("He--", -0.3), ("C+6", 1.0), ("C+6", 20.0))
    sp = Species(name)
    E = Beamlines.R_to_E(sp, R)
    beta_gamma = Beamlines.R_to_beta_gamma(sp, R)
    @test Beamlines.R_to_gamma(sp, R) ≈ E/massof(sp)
    @test Beamlines.R_to_gamma(sp, R) ≈ sqrt(1 + beta_gamma^2)
    @test Beamlines.R_to_v(sp, R) ≈ BT.beta_gamma_to_v(beta_gamma)
    @test Beamlines.R_to_v(sp, R) ≈ Beamlines.E_to_v(sp, E)
    @test Beamlines.R_to_v(sp, R) < BT.C_LIGHT
    @test Beamlines.E_to_R(sp, E) ≈ R
    @test Beamlines.pc_to_R(sp, Beamlines.R_to_pc(sp, R)) ≈ R
  end

  # sincu/sinhcu to machine precision (BigFloat reference) for every argument, including where
  # the series takes over (a series to x^2 below 0.004 was off by 2e-12), and accurate derivatives
  @testset "sincu, sinhcu" begin
    xs = vcat(0.0, 10 .^ range(-12, log10(3.0), length=1001), 0.1 .+ (-3:3) .* eps(0.1), 0.004)
    xs = vcat(xs, -xs)
    ref(f, x) = x == 0 ? big(1.0) : f(big(x))/big(x)
    for T in (Float64, Float32), (f, g) in ((BT.sincu, sin), (BT.sinhcu, sinh))
      @test maximum(x -> abs(f(T(x)) - ref(g, T(x)))/ref(g, T(x)), xs) < eps(T)
    end
    dref(x) = (cos(big(x))*big(x) - sin(big(x)))/big(x)^2
    @test maximum(x -> abs(ForwardDiff.derivative(BT.sincu, x) - dref(x))/abs(dref(x)),
                  filter(!iszero, xs)) < 1e-12
    @test ForwardDiff.derivative(y -> ForwardDiff.derivative(BT.sincu, y), 0.0) ≈ -1/3
  end

  # sincos_quaternion: sin(sqrt(x))/sqrt(x) and cos(sqrt(x)) (sinh, cosh for x < 0) to within eps
  # (BigFloat reference), with the series for |x| < 1 and the functions outside, and derivatives
  # (the previous small-x branch had them 1e-8 off)
  @testset "sincos_quaternion" begin
    sq(x) = x == 0 ? (big(1.0), big(1.0)) : x > 0 ? (sin(sqrt(big(x)))/sqrt(big(x)), cos(sqrt(big(x)))) :
                                                      (sinh(sqrt(-big(x)))/sqrt(-big(x)), cosh(sqrt(-big(x))))
    rel(a, b) = abs(big(a) - b)/abs(b)
    xs = vcat(0.0, 10 .^ range(-12, 0, length=1001), 1 - 1e-12)
    xs = vcat(xs, -xs)
    for T in (Float64, Float32)
      @test maximum(x -> maximum(rel.(BT.sincos_quaternion(T(x)), sq(T(x)))), xs) < 2*eps(T)
      # outside the series: the functions
      @test all(x -> BT.sincos_quaternion(T(x)) === BT._sincos_quaternion(T(x)), (1.0, -1.0, 2.5, -30.0))
    end
    dsq(x, k) = (h = big(1e-30); (sq(x + h)[k] - sq(x - h)[k])/(2h))
    for k in 1:2
      @test maximum(x -> rel(ForwardDiff.derivative(y -> BT.sincos_quaternion(y)[k], x), dsq(x, k)),
                    filter(x -> abs(x) > 1e-300, xs)) < 1e-14
    end
  end

  # Series of the exact bend (asin_s, sincos_s, sin_s, atan2_s) to within eps over their ranges,
  # and the functions themselves outside them and for Duals
  @testset "Bend series" begin
    grid(m) = (g = vcat(0.0, m .* 10 .^ range(-10, 0, length=501), m*(1 - 1e-12)); vcat(g, -g))
    rel(a, b) = abs(big(a) - b)/max(abs(b), big(floatmin(Float64)))
    for T in (Float64, Float32)
      @test maximum(x -> rel(BT.asin_s(T(x)), asin(big(T(x)))), grid(0.05)) < 1.5*eps(T)
      @test maximum(x -> rel(BT.sincos_s(T(x))[1], sin(big(T(x)))), grid(0.1)) < 1.5*eps(T)
      @test maximum(x -> rel(BT.sincos_s(T(x))[2], cos(big(T(x)))), grid(0.1)) < 1.5*eps(T)
      @test maximum(x -> rel(BT.sin_s(T(x)), sin(big(T(x)))), grid(0.1)) < 1.5*eps(T)
      @test maximum(r -> rel(BT.atan2_s(3*T(r), T(3)), atan(big(3*T(r))/3)), grid(0.05)) < 1.5*eps(T)
    end
    @test BT.asin_s(0.3) === asin(0.3) && BT.sincos_s(0.5) === sincos(0.5) && BT.sin_s(0.5) === sin(0.5)
    @test BT.atan2_s(0.1, -1.0) === BT.atan2(0.1, -1.0) && BT.atan2_s(0.1, 1.0) === BT.atan2(0.1, 1.0)
    d = ForwardDiff.Dual(0.01, 1.0)
    @test BT.asin_s(d) === asin(d) && BT.sincos_s(d) === sincos(d) && BT.atan2_s(d, 2d) === BT.atan2(d, 2d)
  end
end
