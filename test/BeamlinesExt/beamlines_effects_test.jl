# Effects applied at the entrance of each slice
# (track!(bunch, bl; effects), BeamTracking.effect_setup/effect_kick!)
struct RecordingEffect
  id::Int
  kick::Float64   # pz kick per metre of slice
  log::Vector{Any}
end

function BeamTracking.effect_setup(e::RecordingEffect, bunch, bl)
  push!(e.log, (:setup, e.id, length(bl.line), bunch.p_over_q_ref))
  return 10*e.id
end

function BeamTracking.effect_kick!(e::RecordingEffect, state, bunch,
                                   ele, sp)
  push!(e.log, (:kick, e.id, state, sp, ele.L, bunch.t_ref))
  bunch.coords.v[:, 6] .+= e.kick * sp.L_slice
  return nothing
end

@testset "Slice effects" begin
  e = Species("electron")
  R0 = -10.0
  v0 = 1e-3 .* randn(8, 6)
  mkline() = [Drift(L=1.0), Quadrupole(Kn1=0.2, L=1.2, n_slices=3),
              Drift(L=0.0), Drift(L=0.5, n_slices=2)]
  function run(effects...; line=mkline())
    b = Bunch(copy(v0); p_over_q_ref=R0, species=e)
    bl = Beamline(line; p_over_q_ref=R0, species_ref=e)
    if isempty(effects)
      track!(b, bl; notify=false)
    else
      track!(b, bl; notify=false, effects=effects)
    end
    return b
  end

  # Slice points, setup and tuple order
  log = Any[]
  b = run(RecordingEffect(1, 0.0, log), RecordingEffect(2, 0.0, log))
  @test log[1] == (:setup, 1, 4, R0)
  @test log[2] == (:setup, 2, 4, R0)
  kicks = log[3:end]
  @test length(kicks) == 2*7
  @test all(k -> k[1] === :kick, kicks)
  @test [k[2] for k in kicks] == repeat([1, 2], 7)
  @test [k[3] for k in kicks] == repeat([10, 20], 7)
  sps = [k[4] for k in kicks[1:2:end]]
  SP = BeamTracking.SlicePoint
  Lq = 1.2/3
  Ld = 0.5/2
  @test sps == [SP(1, 1, 1, 0.0, 1.0, 1.0),
                SP(2, 1, 3, 0*Lq, Lq, 1.2),
                SP(2, 2, 3, 1*Lq, Lq, 1.2),
                SP(2, 3, 3, 2*Lq, Lq, 1.2),
                SP(3, 1, 1, 0.0, 0.0, 0.0),
                SP(4, 1, 2, 0*Ld, Ld, 0.5),
                SP(4, 2, 2, 1*Ld, Ld, 0.5)]
  @test [k[5] for k in kicks[1:2:end]] == [1.0, 1.2, 1.2, 1.2, 0, 0.5, 0.5]
  # The reference time is the one left by the previous slice
  t = [k[6] for k in kicks[1:2:end]]
  @test t[1] == 0 && issorted(t) && t[5] == t[6] && t[7] > t[6]
  @test t[7] < b.t_ref

  # Effects that do nothing, and no effects, are the plain track!
  plain = run()
  @test b.coords.v == plain.coords.v
  @test b.t_ref == plain.t_ref
  b = Bunch(copy(v0); p_over_q_ref=R0, species=e)
  track!(b, Beamline(mkline(); p_over_q_ref=R0, species_ref=e);
         notify=false, effects=())
  @test b.coords.v == plain.coords.v

  # The kick acts before each slice is tracked: same as kicking by
  # hand between separately tracked elements of the slice lengths
  c = 1e-4
  kicked = run(RecordingEffect(1, c, Any[]);
               line=[Drift(L=0.5, n_slices=2), Drift(L=1.0)])
  hand = Bunch(copy(v0); p_over_q_ref=R0, species=e)
  for L in (0.25, 0.25, 1.0)
    hand.coords.v[:, 6] .+= c * L
    track!(hand, Beamline([Drift(L=L)]; p_over_q_ref=R0, species_ref=e);
           notify=false)
  end
  @test maximum(abs.(kicked.coords.v .- hand.coords.v)) <= 1e-15
  @test kicked.t_ref ≈ hand.t_ref rtol=1e-15
  @test maximum(abs.(kicked.coords.v .- plain.coords.v)) > 1e-6

end
