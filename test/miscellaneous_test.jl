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
end
