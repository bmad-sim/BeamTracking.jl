"""
    reference_momentum_shift!(i, coords::Coords, beta_gamma_old, dbeta_gamma, shift_pz)

Shift coordinates due to a change in reference energy `dbeta_gamma`.
"""
@inline function reference_momentum_shift!(i, coords::Coords, beta_gamma_old, dbeta_gamma, ::Val{shift_pz}) where {shift_pz}
  @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)

    beta_gamma_new = beta_gamma_old + dbeta_gamma
    beta_gamma_ratio = beta_gamma_old / beta_gamma_new
    
    v[i,PXI] = vifelse(alive, beta_gamma_ratio * v[i,PXI], v[i,PXI])
    v[i,PYI] = vifelse(alive, beta_gamma_ratio * v[i,PYI], v[i,PYI])
    if shift_pz
      v[i,PZI] = vifelse(alive, fma(beta_gamma_old, v[i,PZI], -dbeta_gamma) / beta_gamma_new, v[i,PZI])
    end
  end
end

# Fallback
fma(a, b, c) = a * b + c



