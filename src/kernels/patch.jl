@makekernel fastgtpsa=true function patch_offset!(i, coords::Coords, tilde_m, dx, dy, dt)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p = 1 + v[i,PZI]
  new_x = v[i,XI] - dx
  new_y = v[i,YI] - dy
  new_z = v[i,ZI] + rel_p/sqrt(rel_p*rel_p+tilde_m*tilde_m)*C_LIGHT*dt
  v[i,XI] = vifelse(alive, new_x, v[i,XI])
  v[i,YI] = vifelse(alive, new_y, v[i,YI])
  v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
end

@makekernel fastgtpsa=true function patch!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, dt, dx, dy, dz, winv, L) 
  v = coords.v
  rel_p = 1 + v[i,PZI]
  ps_02 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  ps_02_0 = zero(ps_02)
  good_momenta = (ps_02 > ps_02_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps_02_1 = one(ps_02)
  ps_0 = sqrt(vifelse(good_momenta, ps_02, ps_02_1))
  # Only apply rotations if needed
  if isnothing(winv)
    patch_offset!(i, coords, tilde_m, dx, dy, dt)
    exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, L)
    new_z = v[i,ZI] - (dz - L) * rel_p / ps_0
    v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
  else
    patch_offset!(i, coords, tilde_m, dx, dy, dt)
    w31 = 2*(winv[QX]*winv[QZ] - winv[QY]*winv[Q0])
    w32 = 2*(winv[QY]*winv[QZ] + winv[QX]*winv[Q0])
    w33 = 1 - 2*(winv[QX]*winv[QX] + winv[QY]*winv[QY])
    s_f = w31*v[i,XI] + w32*v[i,YI] - w33*dz
    rotation!(i, coords, winv, -dz)
    exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, -s_f)
    new_z = v[i,ZI] + ((s_f + L) * rel_p * 
    sqrt((1 + tilde_m*tilde_m)/(rel_p*rel_p + tilde_m*tilde_m)))
    v[i,ZI] = vifelse(alive, new_z, v[i,ZI])
  end
end
