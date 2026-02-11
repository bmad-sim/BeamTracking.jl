"""
    reference_momentum_shift!(i, coords::Coords, P0c_old, dP0c, shift_pz)

Shift coordinates due to a change in reference energy `dE`.
"""
@makekernel fastgtpsa=true function reference_momentum_shift!(i, coords::Coords, P0c_old, dP0c, shift_pz)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  P0c_new = P0c_old + dP0c
  P0c_ratio = P0c_old / P0c_new

  v[i,PXI] = vifelse(alive, P0c_ratio * v[i,PXI], v[i,PXI])
  v[i,PYI] = vifelse(alive, P0c_ratio * v[i,PYI], v[i,PYI])
  if shift_pz
    v[i,PZI] = vifelse(alive, P0c_ratio * v[i,PZI] - dP0c / P0c_new, v[i,PZI])
  end
end
