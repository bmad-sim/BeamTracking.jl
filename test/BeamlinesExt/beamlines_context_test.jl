@testset "Beamlines Context test" begin
    # get sample data
  ele1 = Drift(L = 5)
  ele2 = Drift(L = 6)
  bl1 = Beamline([ele1], E_ref=18e9, species_ref=Species("electron"))
  bl2 = Beamline([ele2], E_ref=18e9, species_ref=Species("electron"))
  b01 = Bunch(collect(transpose(@vars(D1))), nothing, p_over_q_ref=bl1.p_over_q_ref, species=Species("electron"))
  b02 = Bunch(collect(transpose(@vars(D1))), nothing, p_over_q_ref=bl1.p_over_q_ref, species=Species("electron"))
  track!(b01, bl1)
  track!(b02, bl2)
  v1 = GTPSA.jacobian(b01.coords.v)
  v2 = GTPSA.jacobian(b02.coords.v)

  # Now use contexts: the Beamline's own, with values tried temporarily (push!/pop!, withvalues)
  ele = Drift(L=DefExpr(c -> c.L))
  bl = Beamline([ele], E_ref=18e9, species_ref=Species("electron"), context=Context(L = 5))
  newbunch() = Bunch(collect(transpose(@vars(D1))), nothing, p_over_q_ref=bl1.p_over_q_ref,
                     species=Species("electron"))
  jac(bunch) = GTPSA.jacobian(bunch.coords.v)
  b0 = newbunch(); track!(b0, bl)
  @test jac(b0) ≈ v1
  push!(bl.context; L = 6)
  b0 = newbunch(); track!(b0, bl)
  @test jac(b0) ≈ v2
  pop!(bl.context)
  b0 = newbunch(); track!(b0, bl)
  @test jac(b0) ≈ v1
  withvalues(bl.context; L = 6) do
    b0 = newbunch(); track!(b0, bl)
    @test jac(b0) ≈ v2
  end
  @test bl.context.L == 5

  # An element that is not in a Beamline uses DEFAULT_CONTEXT
  elef = Drift(L=DefExpr(c -> c.Lfree))
  withvalues(DEFAULT_CONTEXT; Lfree = 6) do
    b0 = newbunch(); track!(b0, elef)
    @test jac(b0) ≈ v2
  end
end
