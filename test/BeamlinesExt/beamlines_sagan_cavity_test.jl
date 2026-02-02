include("../lattices/sagan_cavity_lat.jl")
import BeamTracking as BT

v1 = [0.01  0.02  0.03  0.04  0.05  0.1]

R0 = BT.E_to_R(species, E0)

@testset "SaganCavity" begin
  ele = sc1
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  @test b1.coords.v ≈ [0.04217873136425864 0.018293913374049727 0.09332500470086044 0.03619832737021669 0.06606835704543228 0.14023826463265046]

  ele = sc2
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  @test b1.coords.v ≈ [0.04194241427584856 0.018231932600906685 0.09388482855169711 0.03646386520181337 0.06713366175549244 0.3270142635510056]

  ele = sc3
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  #@test b1.coords.v ≈ [0.01 0.018382799381485956 0.03 0.03676559876297191 0.04895340854206434 -0.27490000696709455]

  ele = sc4
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  #@test b1.coords.v ≈ [0.01 0.018048002564646365 0.03 0.03609600512929273 0.05054378625481679 0.20674145046402864]

  ele = sc5
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  #@test b1.coords.v ≈ [0.04013979184863006 0.02073739380288904 0.09063188862813992 0.04363580632735068 0.0685983737708154 0.33296212682308773]

  ele = sc6
  b1 = Bunch(deepcopy(v1), species=species, p_over_q_ref = BT.E_to_R(species, ele.E_ref-ele.dE_ref))
  track!(b1, ele)
  #@test b1.coords.v ≈ [0.05411818209039813 0.017031536764595845 0.11975670804339537 0.03399302609195471 0.03855654036685451 -0.2129528915963015]
end







# Test spin
# Test multipoles
# Test zero length
# standing/traveling waves
# dE_ref
# differing n_steps
# differing n_cell
# non-zero starting orbit
