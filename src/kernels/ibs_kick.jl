@makekernel function ibs_damping_and_diffusion!(i, coords::Coords, tilde_m, gamma_0, b_coeff, integrals, diffusion, P, sigma_inv, means, g, w, w_inv, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  h = 1 + g*v[i,XI]

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end

  rel_p = 1 + v[i,PZI]
  beta_gamma = rel_p/tilde_m 
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma/gamma

  pl2 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  dt_ds = h*rel_p/(beta*C_LIGHT*pl)
  moving_forward = (dt_ds > zero(dt_ds))
  coords.state[i] = vifelse(!moving_forward & alive, STATE_LOST, coords.state[i])
  dt_ds = vifelse(moving_forward, dt_ds, one(dt_ds))
  dt = dt_ds*L

  X =   v[i,XI]  - means[XI]
  Y =   v[i,YI]  - means[YI]
  Z =  (v[i,ZI]  - means[ZI] )#*gamma_0
  PX = (v[i,PXI] - means[PXI])#*p0
  PY = (v[i,PYI] - means[PYI])#*p0
  PZ = (v[i,PZI] - means[PZI])#*p0/gamma_0

  kx = X*sigma_inv[1,2] + Y*sigma_inv[3,2] + Z*sigma_inv[5,2] + PX*sigma_inv[2,2] + PY*sigma_inv[4,2] + PZ*sigma_inv[6,2]
  ky = X*sigma_inv[1,4] + Y*sigma_inv[3,4] + Z*sigma_inv[5,4] + PX*sigma_inv[2,4] + PY*sigma_inv[4,4] + PZ*sigma_inv[6,4]
  kz = X*sigma_inv[1,6] + Y*sigma_inv[3,6] + Z*sigma_inv[5,6] + PX*sigma_inv[2,6] + PY*sigma_inv[4,6] + PZ*sigma_inv[6,6]
  kz *= gamma_0
  
  wx = P[1,1]*kx + P[1,2]*ky + P[1,3]*kz
  wy = P[2,1]*kx + P[2,2]*ky + P[2,3]*kz
  wz = P[3,1]*kx + P[3,2]*ky + P[3,3]*kz

  mx = wx*integrals[1]
  my = wy*integrals[2]
  mz = wz*integrals[3]

  I_X = P[1,1]*mx + P[2,1]*my + P[3,1]*mz
  I_Y = P[1,2]*mx + P[2,2]*my + P[3,2]*mz
  I_Z = P[1,3]*mx + P[2,3]*my + P[3,3]*mz

  X_A_X  =   sigma_inv[1,1]*X*X +   sigma_inv[3,3]*Y*Y +   sigma_inv[5,5]*Z*Z
  X_A_X += 2*sigma_inv[1,3]*X*Y + 2*sigma_inv[1,5]*X*Z + 2*sigma_inv[3,5]*Y*Z

  X_B_P  =   X*(sigma_inv[1,2]*PX + sigma_inv[1,4]*PY + sigma_inv[1,6]*PZ)
  X_B_P +=   Y*(sigma_inv[2,3]*PX + sigma_inv[3,4]*PY + sigma_inv[3,6]*PZ) 
  X_B_P +=   Z*(sigma_inv[2,5]*PX + sigma_inv[4,5]*PY + sigma_inv[5,6]*PZ)

  P_C_P  =   sigma_inv[2,2]*PX*PX +   sigma_inv[4,4]*PY*PY +   sigma_inv[6,6]*PZ*PZ
  P_C_P += 2*sigma_inv[2,4]*PX*PY + 2*sigma_inv[2,6]*PX*PZ + 2*sigma_inv[4,6]*PY*PZ

  exp_correction = exp(-X_A_X/2 - X_B_P - P_C_P/2)

  B_X = exp_correction*b_coeff*I_X 
  B_Y = exp_correction*b_coeff*I_Y
  B_Z = exp_correction*b_coeff*I_Z

  b_x = B_X/gamma_0
  b_y = B_Y/gamma_0
  b_z = B_Z

  d_xx, d_yy, d_zz = exp_correction .* diffusion
  #println(d_zz)

  sigma_x = sqrt(d_xx*dt)
  sigma_y = sqrt(d_yy*dt)
  sigma_z = sqrt(d_zz*dt)

  rand_px, rand_py = gaussian_random(sigma_x, sigma_y)
  rand_pz, _       = gaussian_random(sigma_z, zero(sigma_z))

  new_px = v[i,PXI] + b_x*dt + rand_px
  new_py = v[i,PYI] + b_y*dt + rand_py
  new_pz = v[i,PZI] + b_z*dt + rand_pz

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end
#=
@makekernel function ibs_damping_and_diffusion!(i, coords::Coords, tilde_m, gamma_0, p0, b_coeff, integrals, diffusion, P, M_inv, means, g, w, w_inv, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  # I assume we are between elements so there is no vector potential.
  #if mm[1] == 0
  #  ax = -v[i,YI] * kn[1] / 2
  #  ay =  v[i,XI] * kn[1] / 2
  #else
  #  ax = zero(v[i,XI])
  #  ay = ax
  #end

  #px = v[i,PXI] - ax
  #py = v[i,PYI] - ay

  h = 1 + g*v[i,XI]

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end

  rel_p = 1 + v[i,PZI]
  beta_gamma = rel_p/tilde_m 
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma/gamma

  pl2 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  dt_ds = h*rel_p/(beta*C_LIGHT*pl)
  moving_forward = (dt_ds > zero(dt_ds))
  coords.state[i] = vifelse(!moving_forward & alive, STATE_LOST, coords.state[i])
  dt_ds = vifelse(moving_forward, dt_ds, one(dt_ds))
  dt = dt_ds*L

  # Transformation to beam frame and relative to center of beam
  X =   v[i,XI]  - means[XI]
  Y =   v[i,YI]  - means[YI]
  Z =  (v[i,ZI]  - means[ZI] )*gamma_0
  PX = (v[i,PXI] - means[PXI])*p0
  PY = (v[i,PYI] - means[PYI])*p0
  PZ = (v[i,PZI] - means[PZI])*p0/gamma_0

  kx = X*M_inv[1,2] + Y*M_inv[3,2] + Z*M_inv[5,2] + PX*M_inv[2,2] + PY*M_inv[4,2] + PZ*M_inv[6,2]
  ky = X*M_inv[1,4] + Y*M_inv[3,4] + Z*M_inv[5,4] + PX*M_inv[2,4] + PY*M_inv[4,4] + PZ*M_inv[6,4]
  kz = X*M_inv[1,6] + Y*M_inv[3,6] + Z*M_inv[5,6] + PX*M_inv[2,6] + PY*M_inv[4,6] + PZ*M_inv[6,6]
  
  wx = P[1,1]*kx + P[1,2]*ky + P[1,3]*kz
  wy = P[2,1]*kx + P[2,2]*ky + P[2,3]*kz
  wz = P[3,1]*kx + P[3,2]*ky + P[3,3]*kz

  I_X = P[1,1]*wx*integrals[1] + P[2,1]*wy*integrals[2] + P[3,1]*wz*integrals[3]
  I_Y = P[1,2]*wx*integrals[1] + P[2,2]*wy*integrals[2] + P[3,2]*wz*integrals[3]
  I_Z = P[1,3]*wx*integrals[1] + P[2,3]*wy*integrals[2] + P[3,3]*wz*integrals[3]

  B_X = b_coeff*I_X 
  B_Y = b_coeff*I_Y
  B_Z = b_coeff*I_Z

  b_x = B_X/gamma_0/p0
  b_y = B_Y/gamma_0/p0
  b_z = B_Z/p0

  d_xx, d_yy, d_zz = diffusion
  println(b_x, " ", b_y, " ", b_z)

  sigma_x = sqrt(d_xx*dt)
  sigma_y = sqrt(d_yy*dt)
  sigma_z = sqrt(d_zz*dt)

  rand_px, rand_py = gaussian_random(sigma_x, sigma_y)
  rand_pz, _       = gaussian_random(sigma_z, zero(sigma_z))

  new_px = v[i,PXI] + b_x*dt + rand_px
  new_py = v[i,PYI] + b_y*dt + rand_py
  new_pz = v[i,PZI] + b_z*dt + rand_pz

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end
=#