@testset "Branch tracking" begin
  species = Species("electron")
  E1 = 1e9
  dE = 0.5e9
  bl1 = Beamline([Drift(L=1.0), Quadrupole(L=0.5, Kn1=0.2)], E_ref=E1, species_ref=species)
  bl2 = Beamline([Drift(L=2.0), Quadrupole(L=0.5, Kn1=-0.2)], dE_ref=dE)
  branch = Branch([bl1, bl2])

  R1 = branch.beamlines[1].p_over_q_ref
  R2 = branch.beamlines[2].p_over_q_ref
  @test BeamTracking.R_to_E(species, R2) ≈ E1 + dE

  v0 = [1e-3 2e-3 -1e-3 3e-4 1e-2 -2e-3
        -2e-3 -1e-3 2e-3 -5e-4 -1e-2 5e-3]

  # Species and reference energy are set from the first Beamline
  b = Bunch(copy(v0))
  track!(b, branch)
  @test b.species == species
  @test b.p_over_q_ref ≈ R2

  # Expected: track bl1, rescale coordinates to the new reference with the particle
  # momenta unchanged, then track bl2
  b_exp = Bunch(copy(v0), p_over_q_ref=R1, species=species)
  track!(b_exp, Beamline([Drift(L=1.0), Quadrupole(L=0.5, Kn1=0.2)], p_over_q_ref=R1, species_ref=species))
  pc_before = (1 .+ b_exp.v[:,BeamTracking.PZI]) .* R1
  px_before = b_exp.v[:,BeamTracking.PXI] .* R1
  b_exp.v[:,BeamTracking.PXI] .*= R1/R2
  b_exp.v[:,BeamTracking.PYI] .*= R1/R2
  b_exp.v[:,BeamTracking.PZI] .= (1 .+ b_exp.v[:,BeamTracking.PZI]) .* R1 ./ R2 .- 1
  @test (1 .+ b_exp.v[:,BeamTracking.PZI]) .* R2 ≈ pc_before
  @test b_exp.v[:,BeamTracking.PXI] .* R2 ≈ px_before
  b_exp.p_over_q_ref = R2
  track!(b_exp, Beamline([Drift(L=2.0), Quadrupole(L=0.5, Kn1=-0.2)], p_over_q_ref=R2, species_ref=species))
  @test b.v ≈ b_exp.v
  @test b.t_ref ≈ 1.5 / BeamTracking.R_to_v(species, R1) + 2.5 / BeamTracking.R_to_v(species, R2)

  # The step in reference energy between Beamlines never changes the particle energies
  b2 = Bunch(copy(v0))
  track!(b2, branch; ramp_particle_energy_without_rf=true)
  @test b2.v ≈ b.v

  # Same as tracking the Beamlines of the Branch one after the other with a uniform energy
  bl3 = Beamline([Drift(L=2.0)], E_ref=E1, species_ref=species)
  b3 = Bunch(copy(v0))
  b4 = Bunch(copy(v0))
  track!(b3, Branch([bl1, bl3]))
  track!(b4, bl1)
  track!(b4, bl3)
  @test b3.v ≈ b4.v

  # Branch made from LineElements with a change in reference energy
  branch5 = Branch([Marker(E_ref=E1, species_ref=species), Drift(L=1.0), Marker(dE_ref=dE), Drift(L=2.0)])
  b5 = Bunch(copy(v0))
  track!(b5, branch5)
  b6 = Bunch(copy(v0))
  track!(b6, Branch([Beamline([Drift(L=1.0)], E_ref=E1, species_ref=species), Beamline([Drift(L=2.0)], dE_ref=dE)]))
  @test b5.v ≈ b6.v
  @test b5.p_over_q_ref ≈ R2

  # Empty Branch
  b7 = Bunch(copy(v0))
  track!(b7, Branch(Beamline[]))
  @test b7.v == v0
end
