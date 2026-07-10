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

    # Now use contexts
    push!(GLOBAL_CONTEXTS, Context(L = 5))
    c = Context(L = 6)
    ele = Drift(L=DefExpr(c -> c.L))
    bl = Beamline([ele], E_ref=18e9, species_ref=Species("electron"))
    b0 = Bunch(collect(transpose(@vars(D1))), nothing, p_over_q_ref=bl1.p_over_q_ref, species=Species("electron"))
    track!(b0, bl)
    v1c = GTPSA.jacobian(b0.coords.v)
    bl.context = c
    b0 = Bunch(collect(transpose(@vars(D1))), nothing, p_over_q_ref=bl1.p_over_q_ref, species=Species("electron"))
    track!(b0, bl)
    v2c = GTPSA.jacobian(b0.coords.v)
    @test v1 ≈ v1c
    @test v2 ≈ v2c
end