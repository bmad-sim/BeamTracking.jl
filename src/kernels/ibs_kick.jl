@makekernel function ibs_damping_and_diffusion!(i, coords::Coords, backend, tilde_m, gamma_0, ::Val{damping_on}, ::Val{fluctuations_on}, b_coeff, integrals, diffusion_lambdas, diffusion_P, P, sigma_inv, means, g, w, w_inv, L)
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
  Z =   v[i,ZI]  - means[ZI] 
  PX =  v[i,PXI] - means[PXI]
  PY =  v[i,PYI] - means[PYI]
  PZ =  v[i,PZI] - means[PZI]

  kx = sigma_inv[2]*X + sigma_inv[8] *Y + sigma_inv[10]*Z + sigma_inv[7] *PX + sigma_inv[9] *PY + sigma_inv[11]*PZ
  ky = sigma_inv[4]*X + sigma_inv[13]*Y + sigma_inv[17]*Z + sigma_inv[9] *PX + sigma_inv[16]*PY + sigma_inv[18]*PZ
  kz = sigma_inv[6]*X + sigma_inv[15]*Y + sigma_inv[20]*Z + sigma_inv[11]*PX + sigma_inv[18]*PY + sigma_inv[21]*PZ
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

  X_A_X  =    sigma_inv[1]*X*X + sigma_inv[12]*Y*Y + sigma_inv[19]*Z*Z
  X_A_X += 2*(sigma_inv[3]*X*Y + sigma_inv[5]*X*Z  + sigma_inv[14]*Y*Z)

  X_B_P  =   X*(sigma_inv[2]*PX  + sigma_inv[4]*PY  + sigma_inv[6]*PZ)
  X_B_P +=   Y*(sigma_inv[8]*PX  + sigma_inv[13]*PY + sigma_inv[15]*PZ) 
  X_B_P +=   Z*(sigma_inv[10]*PX + sigma_inv[17]*PY + sigma_inv[20]*PZ)

  P_C_P  =    sigma_inv[7]*PX*PX + sigma_inv[16]*PY*PY + sigma_inv[21]*PZ*PZ
  P_C_P += 2*(sigma_inv[9]*PX*PY + sigma_inv[11]*PX*PZ + sigma_inv[18]*PY*PZ)

  exp_correction = 2*exp(-X_A_X/2 - X_B_P - P_C_P/2)

  B_X = exp_correction*b_coeff*I_X 
  B_Y = exp_correction*b_coeff*I_Y
  B_Z = exp_correction*b_coeff*I_Z

  b_x = B_X/gamma_0
  b_y = B_Y/gamma_0
  b_z = B_Z

  new_px = v[i,PXI] + b_x*dt
  new_py = v[i,PYI] + b_y*dt
  new_pz = v[i,PZI] + b_z*dt

  v[i,PXI] = vifelse(alive && damping_on, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive && damping_on, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive && damping_on, new_pz, v[i,PZI])

  # Rotate in with P

  new_px = diffusion_P[1,1]*v[i,PXI] + diffusion_P[1,2]*v[i,PYI] + diffusion_P[1,3]*v[i,PZI]
  new_py = diffusion_P[2,1]*v[i,PXI] + diffusion_P[2,2]*v[i,PYI] + diffusion_P[2,3]*v[i,PZI]
  new_pz = diffusion_P[3,1]*v[i,PXI] + diffusion_P[3,2]*v[i,PYI] + diffusion_P[3,3]*v[i,PZI]

  v[i,PXI] = vifelse(alive && fluctuations_on, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive && fluctuations_on, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive && fluctuations_on, new_pz, v[i,PZI])

  d_xx, d_yy, d_zz = exp_correction .* diffusion_lambdas

  sigma_x = sqrt(d_xx*dt)
  sigma_y = sqrt(d_yy*dt)
  sigma_z = sqrt(d_zz*dt)

  rand_px, rand_py = gaussian_random(backend, sigma_x, sigma_y)
  rand_pz, _       = gaussian_random(backend, sigma_z, zero(sigma_z))

  new_px = v[i,PXI] + rand_px
  new_py = v[i,PYI] + rand_py
  new_pz = v[i,PZI] + rand_pz

  v[i,PXI] = vifelse(alive && fluctuations_on, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive && fluctuations_on, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive && fluctuations_on, new_pz, v[i,PZI])

  # Rotate out with P'

  new_px = diffusion_P[1,1]*v[i,PXI] + diffusion_P[2,1]*v[i,PYI] + diffusion_P[3,1]*v[i,PZI]
  new_py = diffusion_P[1,2]*v[i,PXI] + diffusion_P[2,2]*v[i,PYI] + diffusion_P[3,2]*v[i,PZI]
  new_pz = diffusion_P[1,3]*v[i,PXI] + diffusion_P[2,3]*v[i,PYI] + diffusion_P[3,3]*v[i,PZI]

  v[i,PXI] = vifelse(alive && fluctuations_on, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive && fluctuations_on, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive && fluctuations_on, new_pz, v[i,PZI])
end