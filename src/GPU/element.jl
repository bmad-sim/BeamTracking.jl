module TrackingLattice
  using ..BeamTracking
  using ..BeamTracking: MatrixKick

  export Element

  const Element = Union{MatrixKick.Drift{Float64}, MatrixKick.Quadrupole{Float64}}

end # TrackingLattice
