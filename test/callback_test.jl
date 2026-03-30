@testset "Callbacks" begin
    ring = include("lattices/esr.jl")
    foreach(x->x.tracking_method=Yoshida(num_steps=10), ring.line)
    n_thickeles = count(x->x.L != 0, ring.line)

    # One before everything
    s_pos = zeros(1+10*n_thickeles)
    cur_idx = 1
    function s_in_ele(coords, ds_step, g)
        global cur_idx += 1
        global s_pos[cur_idx] = s_pos[cur_idx-1] + ds_step
    end
    b0 = Bunch(v=zeros(1,6), callbacks=(s_in_ele,))
    track!(b0, ring)
    @test s_pos[1] == 0
    @test s_pos[end] ≈ ring.line[end].s_downstream
end