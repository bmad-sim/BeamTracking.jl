@testset "Callbacks" begin
    ring = include("lattices/esr.jl")
    foreach(x->x.tracking_method=Symplectic(n_steps=10), ring.line)
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

    # Straight misalignment
    n_steps = 100
    d = Drift(L=1.2, tracking_method=Symplectic(n_steps=n_steps))
    dx = Drift(L=1.2, tracking_method=Symplectic(n_steps=n_steps), x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6,)
    b  = Beamline([d], species_ref=Species("electron"), E_ref=18e9)
    bx = Beamline([dx], species_ref=Species("electron"), E_ref=18e9)

    orbits = zeros(n_steps, 6)
    quats = zeros(n_steps, 4)
    s = zeros(n_steps)
    t = zeros(n_steps)
    orbitsx = zeros(n_steps, 6)
    quatsx = zeros(n_steps, 4)
    sx = zeros(n_steps)
    tx = zeros(n_steps)
    function savestuff!(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        i_step = round(Int, cur_s/ds_step)
        orbits[i_step, :] = coords.v[1,:]
        quats[i_step, :] = coords.q[1,:]
        s[i_step] = cur_s
        t[i_step] = cur_t_ref
    end
    function savestuffx!(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        i_step = round(Int, cur_s/ds_step)
        transforms_out!(i, coords, cur_s, cur_t_ref)
        orbitsx[i_step, :] = coords.v[1,:]
        quatsx[i_step, :] = coords.q[1,:]
        sx[i_step] = cur_s
        tx[i_step] = cur_t_ref
        transforms_in!(i, coords, cur_s, cur_t_ref)
    end
    b0 = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuff!,))
    b0x = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0, b)
    track!(b0x, bx)
    @test all(orbitsx .≈ orbits)
    @test all(isapprox.(quats, quatsx, atol=1e-16))
    @test all(s .≈ sx)
    @test all(diff(t) .≈ diff(t)[1])
    @test all(diff(tx) .≈ diff(tx)[1])
    @test all(tx .≈ t)

    # Bend misalignment
    n_steps = 100
    d = SBend(L=1.2, g_ref=1e-3, tracking_method=Symplectic(n_steps=n_steps))
    dx = SBend(L=1.2, g_ref=1e-3, tracking_method=Symplectic(n_steps=n_steps), x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6,)
    b  = Beamline([d], species_ref=Species("electron"), E_ref=18e9)
    bx = Beamline([dx], species_ref=Species("electron"), E_ref=18e9)
    b0 = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuff!,))
    b0x = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0, b)
    track!(b0x, bx)
    @test all(orbitsx .≈ orbits)
    @test all(isapprox.(quats, quatsx, atol=1e-13))
    @test all(s .≈ sx)
    @test all(diff(t) .≈ diff(t)[1])
    @test all(diff(tx) .≈ diff(tx)[1])
    @test all(tx .≈ t)


    # Implicit at low energy
    function sol(x, y, s, t, p)
        Ksol= p[1]
        potential = (0.0, -Ksol*y/2, Ksol*x/2, 0.0)
        jac = (0.0,    0.0,     0.0, 0.0,
            0.0,    -Ksol/2, 0.0, 0.0,
            Ksol/2, 0.0,     0.0, 0.0,
            0.0,    0.0,     0.0, 0.0)
        return potential, jac
    end
    solele = Solenoid(four_potential=sol, four_potential_params=(0.3,), four_potential_normalized=true, L=1.2, tracking_method=Symplectic(n_steps=n_steps))
    bx  = Beamline([solele], species_ref=Species("electron"), E_ref=2*massof(Species("electron")))
    b0x = Bunch(v=[0. 0. 0. 0. 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0x, bx)
    @test all(diff(orbitsx[:,5]) .≈ diff(orbitsx[:,5])[1])
    @test all(orbitsx[:,6] .≈ 0.06)

    # Implicit low energy misaligned drift
    function idrift(x, y, s, t, p)
        potential = (0.0, 0.0, 0.0, 0.0)
        jac = (0.0,    0.0,     0.0, 0.0,
            0.0,    0.0, 0.0, 0.0,
            0.0, 0.0,     0.0, 0.0,
            0.0,    0.0,     0.0, 0.0)
        return potential, jac
    end
    d = Drift(L=1.2, tracking_method=Symplectic(n_steps=n_steps), x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6,)
    dx = Drift(L=1.2, four_potential=idrift, tracking_method=Symplectic(n_steps=n_steps), x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6,)
    b  = Beamline([d], species_ref=Species("electron"), E_ref=2*massof(Species("electron")))
    bx = Beamline([dx], species_ref=Species("electron"), E_ref=2*massof(Species("electron")))
    b0 = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0, b)
    orbits .= orbitsx
    quats .= quatsx
    s .= sx
    t .= tx
    b0x = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0x, bx)
    @test all(orbitsx .≈ orbits)
    @test all(isapprox.(quats, quatsx, atol=1e-13))
    @test all(s .≈ sx)
    @test all(diff(t) .≈ diff(t)[1])
    @test all(diff(tx) .≈ diff(tx)[1])
    @test all(tx .≈ t)

    # Implicit crazy (to check nonzero phi works ok)
    # Just to check it works
    function crazy(x, y, s, t, p)
        a0, c0 = p
        potential = (a0*C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), a0*cosh(x)*sin(y)*cos(s)*sinh(t), a0*sinh(x)*cosh(y)*sin(s)*cos(t), a0*c0*cos(x)*sinh(y)*cosh(s)*sin(t))
        jac =  (a0*C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -a0*C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), a0*C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), a0*C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
                a0*cos(s)*sin(y)*sinh(t)*sinh(x),  a0*cos(s)*cos(y)*cosh(x)*sinh(t), -a0*cosh(x)*sin(s)*sin(y)*sinh(t),  a0*cos(s)*cosh(t)*cosh(x)*sin(y),
                a0*cos(t)*cosh(x)*cosh(y)*sin(s),  a0*cos(t)*sin(s)*sinh(x)*sinh(y),  a0*cos(s)*cos(t)*cosh(y)*sinh(x), -a0*cosh(y)*sin(s)*sin(t)*sinh(x),
                -a0*c0*cosh(s)*sin(t)*sin(x)*sinh(y), a0*c0*cos(x)*cosh(s)*cosh(y)*sin(t), a0*c0*cos(x)*sin(t)*sinh(s)*sinh(y), a0*c0*cos(t)*cos(x)*cosh(s)*sinh(y))
        return potential, jac
    end
    crazyelex = Marker(four_potential=crazy, four_potential_params=(1.0, 1.0), four_potential_normalized=true, L=1.2, x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6, tracking_method=Symplectic(n_steps=n_steps))
    bx = Beamline([crazyelex], species_ref=Species("electron"), E_ref=1.1*massof(Species("electron")))
    b0x = Bunch(v=[0.01 0.02 0.03 0.04 0.05 0.06], q=[1.0 0.0 0.0 0.0], callbacks=(savestuffx!,))
    track!(b0x, bx)

    # CUDA test
    #=
    n_steps = 100
    function crazy(x, y, s, t, p)
        a0, c0 = p
        potential = (a0*Beamlines.C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), a0*cosh(x)*sin(y)*cos(s)*sinh(t), a0*sinh(x)*cosh(y)*sin(s)*cos(t), a0*c0*cos(x)*sinh(y)*cosh(s)*sin(t))
        jac =  (a0*Beamlines.C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -a0*Beamlines.C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), a0*Beamlines.C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), a0*Beamlines.C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
                a0*cos(s)*sin(y)*sinh(t)*sinh(x),  a0*cos(s)*cos(y)*cosh(x)*sinh(t), -a0*cosh(x)*sin(s)*sin(y)*sinh(t),  a0*cos(s)*cosh(t)*cosh(x)*sin(y),
                a0*cos(t)*cosh(x)*cosh(y)*sin(s),  a0*cos(t)*sin(s)*sinh(x)*sinh(y),  a0*cos(s)*cos(t)*cosh(y)*sinh(x), -a0*cosh(y)*sin(s)*sin(t)*sinh(x),
                -a0*c0*cosh(s)*sin(t)*sin(x)*sinh(y), a0*c0*cos(x)*cosh(s)*cosh(y)*sin(t), a0*c0*cos(x)*sin(t)*sinh(s)*sinh(y), a0*c0*cos(t)*cos(x)*cosh(s)*sinh(y))
        return potential, jac
    end
    function cusayhi(i, coords, cur_s, cur_t_ref, ds_step, g, transforms_out!, transforms_in!)
        transforms_out!(i, coords, cur_s, cur_t_ref)
        @cuprintln("hello! orbit:", coords.v[i,1], " ", coords.v[i,2], " ", coords.v[i,3], " ", coords.v[i,4], " ", coords.v[i,5], " ", coords.v[i,6], " ")
        @cuprintln("hello! quat:", coords.q[i,1], " ", coords.q[i,2], " ", coords.q[i,3], " ", coords.q[i,4])
        transforms_in!(i, coords, cur_s, cur_t_ref)
    end
    crazyelex = Marker(four_potential=crazy, four_potential_params=(1.0, 1.0), four_potential_normalized=true, L=1.2, x_offset=0.1, y_offset=0.2, z_offset=0.3, x_rot=0.4, y_rot=0.5, tilt=0.6, tracking_method=Symplectic(n_steps=n_steps))
    bx = Beamline([crazyelex], species_ref=Species("electron"), E_ref=1.1*Beamlines.massof(Species("electron")))
    b0x = Bunch(v=CuArray([0.01 0.02 0.03 0.04 0.05 0.06]), q=CuArray([1.0 0.0 0.0 0.0]), callbacks=(cusayhi,))
    track!(b0x, bx)
    =#
end