"""
    reference_momentum_shift!(i, coords::Coords, beta_gamma_old, dbeta_gamma, shift_pz)

Shift coordinates due to a change in reference energy `dbeta_gamma`.
"""
@inline function reference_momentum_shift!(i, coords::Coords, beta_gamma_old, dbeta_gamma, ::Val{shift_pz}) where {shift_pz}
  # No shift (e.g. a prepared chain with a time dependent reference energy that did not change at
  # this element): (beta_gamma_old*pz - 0)/beta_gamma_new would round pz
  all(iszero(dbeta_gamma)) && return nothing
  @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)

    beta_gamma_new = beta_gamma_old + dbeta_gamma
    beta_gamma_ratio = beta_gamma_old / beta_gamma_new
    
    v[i,PXI] = vifelse(alive, beta_gamma_ratio * v[i,PXI], v[i,PXI])
    v[i,PYI] = vifelse(alive, beta_gamma_ratio * v[i,PYI], v[i,PYI])
    if shift_pz
      v[i,PZI] = vifelse(alive, muladd(beta_gamma_old, v[i,PZI], -dbeta_gamma) / beta_gamma_new, v[i,PZI])
    end
  end
end

@inline function reference_time_shift!(i, coords::Coords, tilde_m, c_dt)
  @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)
    rel_p = 1 + v[i,PZI]
    beta = @FastGTPSA rel_p / sqrt(rel_p*rel_p + tilde_m*tilde_m)
    v[i,ZI] = vifelse(alive, muladd(beta, c_dt, v[i,ZI]), v[i,ZI])
  end
end



