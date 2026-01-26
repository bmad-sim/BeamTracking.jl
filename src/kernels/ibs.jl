@makekernel fastgtpsa=true function ibs_damping_and_diffusion!(i, coords::Coords, gamma_0, p0, N, k, L_C, m, det_M, integrals, P, M_inv, means, g, tilt_ref, mm, kn, ks, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  w = rot_quaternion(0, 0, tilt_ref)
  rotation!(i, coords, w, 0)

  if length(mm) > 0 && mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
  else
    ax = zero(v[i,XI])
    ay = ax
  end

  h = 1 + g*v[i,XI]

  w_inv = inv_rot_quaternion(0, 0, tilt_ref)
  rotation!(i, coords, w, 0)

  a = quat_rotate((ax, ay, 0), w_inv)
  px = v[i,PXI] - a[1]
  py = v[i,PYI] - a[2]
  rel_p = 1 + v[i,PZI]

  pl2 = rel_p*rel_p - px*px - py*py
  pl2_0 = zero(pl2)
  good_momenta = (pl2 > pl2_0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pl2_1 = one(pl2)
  pl = sqrt(vifelse(good_momenta, pl2, pl2_1)) 

  dt_ds = h * rel_p / pl
  dt = dt_ds*L

  # Transformation to beam frame and relative to center of beam
  X = v[i,XI] - means[XI]
  Y = v[i,YI] - means[YI]
  Z = (v[i,ZI] - means[ZI])*gamma_0
  PX = (px - means[PXI])*p0
  PY = (py - means[PYI])*p0
  PZ = (v[i,PZI] - means[PZI])*p0*gamma_0

  cx = (P[1,1], P[2,1], P[3,1])
  cy = (P[1,2], P[2,2], P[3,2])
  cz = (P[1,3], P[2,3], P[3,3])

  kx = X*M_inv[1,2] + Y*M_inv[3,2] + Z*M_inv[5,2] + PX*M_inv[2,2] + PY*M_inv[4,2] + PZ*M_inv[6,2]
  ky = X*M_inv[1,4] + Y*M_inv[3,4] + Z*M_inv[5,4] + PX*M_inv[2,4] + PY*M_inv[4,4] + PZ*M_inv[6,4]
  kz = X*M_inv[1,6] + Y*M_inv[3,6] + Z*M_inv[5,6] + PX*M_inv[2,6] + PY*M_inv[4,6] + PZ*M_inv[6,6]
  
  wx = P[1,1] * kx + P[1,2] * ky + P[1,3] * kz
  wy = P[2,1] * kx + P[2,2] * ky + P[2,3] * kz
  wz = P[3,1] * kx + P[3,2] * ky + P[3,3] * kz
  w = (wx, wy, wz)

  I_X = cx[1]*w[1]*integrals[1] + cx[2]*w[2]*integrals[2] + cx[3]*w[3]*integrals[3]
  I_Y = cy[1]*w[1]*integrals[1] + cy[2]*w[2]*integrals[2] + cy[3]*w[3]*integrals[3]
  I_Z = cz[1]*w[1]*integrals[1] + cz[2]*w[2]*integrals[2] + cz[3]*w[3]*integrals[3]

  b_coeff = -N*k*k*L_C/(pi*pi*m*sqrt(det_M))
  B_X = b_coeff*I_X 
  B_Y = b_coeff*I_Y
  B_Z = b_coeff*I_Z

  b_x = B_X/gamma_0/p0
  b_y = B_Y/gamma_0/p0
  b_z = B_Z/p0

  I_XX = (1-cx[1]*cx[1])*integrals[1] + (1-cx[2]*cx[2])*integrals[2] + (1-cx[3]*cx[3])*integrals[3]
  I_YY = (1-cy[1]*cy[1])*integrals[1] + (1-cy[2]*cy[2])*integrals[2] + (1-cy[3]*cy[3])*integrals[3]
  I_ZZ = (1-cz[1]*cz[1])*integrals[1] + (1-cz[2]*cz[2])*integrals[2] + (1-cz[3]*cz[3])*integrals[3]

  d_coeff = -b_coeff/2
  D_XX = d_coeff*I_XX
  D_YY = d_coeff*I_YY
  D_ZZ = d_coeff*I_ZZ

  d_xx = D_XX/gamma_0/p0/p0
  d_yy = D_YY/gamma_0/p0/p0
  d_zz = D_ZZ*gamma_0/p0/p0
  #println(d_zz)

  sigma_x2 = d_xx*dt
  sigma_1 = one(sigma_x2)
  sigma_x = sqrt(vifelse(alive, sigma_x2, sigma_1)) 
  sigma_y2 = d_yy*dt
  sigma_y = sqrt(vifelse(alive, sigma_y2, sigma_1)) 
  sigma_z2 = d_zz*dt
  sigma_z = sqrt(vifelse(alive, sigma_z2, sigma_1))

  new_px = v[i,PXI] + b_x*dt + gaussian_random(sigma_x)
  new_py = v[i,PYI] + b_y*dt + gaussian_random(sigma_y)
  new_pz = v[i,PZI] + b_z*dt + gaussian_random(sigma_z)

  println(new_px - v[i,PXI])
  println(new_py - v[i,PYI])
  println(new_pz - v[i,PZI])

  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
end