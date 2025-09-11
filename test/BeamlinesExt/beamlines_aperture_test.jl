include("../lattices/aperture_lat.jl")

v1 = [ 
       0.0  0.0  6.0  0.0  0.0  0.0
       1.9  0.8  4.0  0.0  0.0  0.0
       2.1 -0.8  4.0  0.0  0.0  0.0
       0.0  0.0  4.0  0.0  0.0  0.0
       0.0  0.0  2.0  0.0  0.0  0.0
       1.5  0.0 -6.0  0.0  0.0  0.0
     ]

v2 = [ 
       1.9  0.0  4.9  0.0  0.0  0.0
       1.8  0.9  4.2  0.0  0.0  0.0
     ]

@testset "Aperture" begin
  @testset "No patch" begin
    bunch_error = Bunch(deepcopy(v1), species=Species("electron"))
    @test_throws ErrorException track!(bunch_error, b_error)
  end

  @testset "Rectangular" begin
    b1 = Bunch(deepcopy(v1), species=Species("electron"))
    track!(b1, b1r)
    @test b1.coords.state == [STATE_ALIVE, STATE_ALIVE, STATE_ALIVE, STATE_ALIVE, STATE_ALIVE, STATE_ALIVE]

    b2 = Bunch(deepcopy(v1), species=Species("electron"))
    track!(b2, b2r)
    @test b2.coords.state == [STATE_LOST, STATE_LOST, STATE_ALIVE, STATE_ALIVE, STATE_LOST, STATE_LOST]

    b3 = Bunch(deepcopy(v1), species=Species("electron"))
    track!(b3, b3r)
    @test b3.coords.state == [STATE_LOST, STATE_LOST, STATE_LOST, STATE_ALIVE, STATE_LOST, STATE_LOST]

    b4 = Bunch(deepcopy(v1), species=Species("electron"))
    track!(b4, b4r)
    @test b4.coords.state == [STATE_LOST, STATE_ALIVE, STATE_LOST, STATE_LOST, STATE_LOST, STATE_ALIVE]
  end

  @testset "Elliptical" begin
    b1 = Bunch(deepcopy(v2), species=Species("electron"))
    track!(b1, b1e)
    @test b1.coords.state == [STATE_ALIVE, STATE_ALIVE]

    b2 = Bunch(deepcopy(v2), species=Species("electron"))
    track!(b2, b2e)
    @test b2.coords.state == [STATE_LOST, STATE_LOST]

    b3 = Bunch(deepcopy(v2), species=Species("electron"))
    track!(b3, b3e)
    @test b3.coords.state == [STATE_LOST, STATE_ALIVE]

    b4 = Bunch(deepcopy(v2), species=Species("electron"))
    @test_throws ErrorException track!(b4, b4e)
  end
end