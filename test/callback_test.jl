@testset "Callbacks" begin
    ring = include("lattices/esr.jl")
    foreach(x->x.tracking_method=Yoshida(n_steps=10), ring.line)
    n_thickeles = count(x->x.L != 0, ring.line)
    n_thineles = count(x->x.L == 0, ring.line)

    # One before everything
    s_pos = zeros(1 + 10*n_thickeles + n_thineles)
    cur_idx = 1
    function s_in_ele(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        cur_idx += 1
        s_pos[cur_idx] = s_pos[cur_idx-1] + ds_step
    end
    b0 = Bunch(v=zeros(1,6), callbacks=(s_in_ele,))
    track!(b0, ring)
    @test s_pos[1] == 0
    @test s_pos[end] ≈ ring.line[end].s_downstream

    d = Drift(L=1.2)
    dx = Drift(L=1.2, x_offset=1, y_offset=2, z_offset=3, x_rot=4, y_rot=5, tilt=6)
    b  = Beamline([d], species_ref=Species("electron"), E_ref=18e9)
    bx = Beamline([dx], species_ref=Species("electron"), E_ref=18e9)

    n_steps = 100
    orbits = zeros(n_steps, 6)
    quats = zeros(n_steps, 4)
    s = zeros(n_steps)
    orbitsx = zeros(n_steps, 6)
    quatsx = zeros(n_steps, 4)
    sx = zeros(n_steps)
    function savestuff!(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        i_step = Int(cur_s/ds_step)
        orbits[i_step, :] = coords.v[1,:]
        quats[i_step, :] = coords.q[1,:]
        s[i_step] = cur_s
    end
    function savestuffx!(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        i_step = Int(cur_s/ds_step)
        @show transforms_out!
        transforms_out!(i, coords, cur_s, cur_t_ref)
        orbitsx[i_step, :] = coords.v[1,:]
        quatsx[i_step, :] = coords.q[1,:]
        sx[i_step] = cur_s
        transforms_in!(i, coords, cur_s, cur_t_ref)
    end
    b0 = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuff!,))
    b0x = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0, b)
    track!(b0x, bx)
end