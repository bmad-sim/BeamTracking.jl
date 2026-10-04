using LinearAlgebra: norm

@testset "RungeKutta hard edge fringe" begin
  species = Species("electron")
  p_over_q_ref = BeamTracking.E_to_R(species, 5e9)
  p0 = [0.01 -0.02 0.015 0.01 0.003 0.03]

  function track_ele(ele)
    bl = Beamline([ele], p_over_q_ref=p_over_q_ref, species_ref=species)
    b = Bunch(copy(p0), p_over_q_ref=p_over_q_ref, species=species)
    track!(b, bl)
    return vec(b.coords.v)
  end

  function apply_fringe!(v, args)
    coords = Coords(fill(STATE_ALIVE, 1), v, nothing, nothing, ())
    BeamTracking.launch!(coords, make_kernel_call(BeamTracking.fringe!, args); use_KA=false)
    return v
  end

  @testset "Constructor" begin
    @test RungeKutta().fringe_at == Fringe.BothEnds
    @test !RungeKutta().multipole_fringe_on
    rk = RungeKutta(n_steps=5, fringe_at=Fringe.ExitEnd, multipole_fringe_on=true)
    @test rk.n_steps == 5 && rk.fringe_at == Fringe.ExitEnd && rk.multipole_fringe_on
    @test RungeKutta(ds_step=0.1, fringe_at=Fringe.NoEnd).fringe_at == Fringe.NoEnd
    @test RungeKutta(0.1, -1) == RungeKutta(ds_step=0.1)
  end

  @testset "Ends are applied around the body" begin
    # Quadrupole: the straight fringe is the multipole fringe kick only
    k1 = 0.8
    quad(fa) = LineElement(L=0.5, Kn1=k1, tracking_method=RungeKutta(n_steps=20, fringe_at=fa))
    rk(fa) = track_ele(quad(fa))
    tilde_m, _, _ = BeamTracking.drift_params(species, p_over_q_ref)
    a = BeamTracking.gyromagnetic_anomaly(species)
    edge = (a, tilde_m, nothing, nothing, nothing, nothing, SA[2], SA[k1], SA[0.0])

    # Entrance fringe, then the body without fringes, then the exit fringe
    v = apply_fringe!(copy(p0), (edge..., 1))
    b = Bunch(v, p_over_q_ref=p_over_q_ref, species=species)
    track!(b, Beamline([quad(Fringe.NoEnd)], p_over_q_ref=p_over_q_ref, species_ref=species))
    @test vec(b.coords.v) ≈ rk(Fringe.EntranceEnd) rtol=1e-14
    @test vec(apply_fringe!(copy(b.coords.v), (edge..., -1))) ≈ rk(Fringe.BothEnds) rtol=1e-14
    @test vec(apply_fringe!(reshape(rk(Fringe.NoEnd), 1, 6), (edge..., -1))) ≈ rk(Fringe.ExitEnd) rtol=1e-14
    @test norm(rk(Fringe.BothEnds) - rk(Fringe.NoEnd)) > 1e-8
  end

  # The fringe kicks are the same as for the symplectic methods so, where the body models agree,
  # RungeKutta with enough steps agrees with symplectic tracking.
  @testset "Agreement with symplectic tracking" begin
    all_ends = (Fringe.BothEnds, Fringe.EntranceEnd, Fringe.ExitEnd)
    bend = (L=1.0, g=0.1, Kn0=0.12, edge1_int=0.02, edge2_int=0.01)
    cases = (
      ((L=0.5, Kn1=0.8), MatrixKick, (;), all_ends),
      ((L=0.5, Kn0=0.02, Ks0=0.01), DriftKick, (;), all_ends),
      ((L=0.3, Kn2=20.0, Ks3=-300.0), DriftKick, (multipole_fringe_on=true,), all_ends),
      # Without the fringe, the RungeKutta solenoid field ends abruptly while the symplectic
      # methods keep the linear edge focusing from the canonical momenta.
      ((L=0.5, Ksol=0.3, Kn1=0.4, Kn2=10.0), SolenoidKick, (multipole_fringe_on=true,), (Fringe.BothEnds,)),
      ((; bend..., e1=0.1, e2=-0.05), BendKick, (;), (Fringe.BothEnds,)),
      ((; bend..., e1=0.1), BendKick, (;), (Fringe.EntranceEnd,)),
      ((; bend..., e2=-0.05), BendKick, (;), (Fringe.ExitEnd,)),
    )
    for (fields, TM, kw, ends) in cases, fa in ends
      v_rk = track_ele(LineElement(; fields..., tracking_method=RungeKutta(; n_steps=400, fringe_at=fa, kw...)))
      v_sy = track_ele(LineElement(; fields..., tracking_method=TM(; order=8, n_steps=20, fringe_at=fa, kw...)))
      @test maximum(abs, v_rk - v_sy) < 1e-12
    end
  end

  @testset "multipole_fringe_on" begin
    rk(on, fa) = track_ele(LineElement(L=0.3, Kn2=20.0,
                                       tracking_method=RungeKutta(n_steps=20, fringe_at=fa, multipole_fringe_on=on)))
    @test rk(false, Fringe.BothEnds) == rk(false, Fringe.NoEnd)
    @test norm(rk(true, Fringe.BothEnds) - rk(true, Fringe.NoEnd)) > 1e-10
  end

  @testset "Bend edge angles require the fringe at that end" begin
    bend(fa; e1=0.0, e2=0.0) = LineElement(L=1.0, g=0.1, Kn0=0.1, e1=e1, e2=e2,
                                           tracking_method=RungeKutta(n_steps=10, fringe_at=fa))
    @test_throws ErrorException track_ele(bend(Fringe.ExitEnd, e1=0.01))
    @test_throws ErrorException track_ele(bend(Fringe.NoEnd, e1=0.01))
    @test_throws ErrorException track_ele(bend(Fringe.EntranceEnd, e2=0.01))
    @test_throws ErrorException track_ele(bend(Fringe.NoEnd, e2=0.01))
    @test all(isfinite, track_ele(bend(Fringe.EntranceEnd, e1=0.01)))
    @test all(isfinite, track_ele(bend(Fringe.ExitEnd, e2=0.01)))
    @test all(isfinite, track_ele(bend(Fringe.BothEnds, e1=0.01, e2=0.01)))
  end
end
