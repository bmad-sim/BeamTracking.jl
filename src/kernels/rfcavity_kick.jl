@makekernel fastgtpsa=true function cavity!(i, coords::Coords, q, mc2, radiation_damping, beta_0, gamsqr_0, tilde_m, E_ref, p0c, a, omega, t0, E0_over_Rref, mm, kn, ks, L)
  multipoles = (length(mm) > 0)
  sol = (multipoles && mm[1] == 0)
  if sol
    exact_solenoid!(i, coords, kn[1], beta_0, gamsqr_0, tilde_m, L/2)
  else
    exact_drift!(   i, coords, beta_0, gamsqr_0, tilde_m, L/2)
  end
  #t0 = t0 + (L/2)/(beta_0*C_LIGHT)

  if multipoles
    if radiation_damping
      deterministic_radiation!(   i, coords, q, mc2, E_ref, 0, mm, kn, ks, L/2)
    end
    multipole_kick!(i, coords, mm, kn * L/2, ks * L/2, -1)
  end

  if isnothing(coords.q)
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, t0, E0_over_Rref, L)
  else
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, t0, E0_over_Rref, L/2)
    rotate_spin_cavity!(          i, coords, a, tilde_m, omega, t0, E0_over_Rref, mm, kn, ks, L)
    cavity_kick!(                 i, coords, beta_0, tilde_m, E_ref, p0c, omega, t0, E0_over_Rref, L/2)
  end

  if multipoles
    multipole_kick!(i, coords, mm, kn * L/2, ks * L/2, -1)
    if radiation_damping
      deterministic_radiation!(   i, coords, q, mc2, E_ref, 0, mm, kn, ks, L/2)
    end
  end

  if sol
    exact_solenoid!(i, coords, kn[1], beta_0, gamsqr_0, tilde_m, L/2)
  else
    exact_drift!(   i, coords, beta_0, gamsqr_0, tilde_m, L/2)
  end
end


@makekernel fastgtpsa=true function bmad_to_mad!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c)
  v = coords.v

  rel_p = 1 + v[i,PZI]
  beta_gamma = rel_p/tilde_m
  gamma = sqrt(1 + beta_gamma*beta_gamma)
  beta = beta_gamma/gamma
  tau = v[i,ZI]/beta

  gamma_0_inv = tilde_m*beta_0
  E = E_ref*gamma*gamma_0_inv

  v[i,ZI]  = tau
  v[i,PZI] = E/p0c - 1/beta_0
end


@makekernel fastgtpsa=true function mad_to_bmad!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c)
  v = coords.v

  E = E_ref + p0c*v[i,PZI]
  gamma_0_inv = tilde_m*beta_0
  gamma = E/E_ref/gamma_0_inv
  beta = sqrt(1-1/(gamma*gamma))
  z = v[i,ZI]*beta
  
  pc = beta*E

  v[i,ZI]  =  z
  v[i,PZI] = (pc-p0c)/p0c
end


@makekernel fastgtpsa=true function cavity_kick!(i, coords::Coords, beta_0, tilde_m, E_ref, p0c, omega, t0, E0_over_Rref, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  bmad_to_mad!(i, coords, beta_0, tilde_m, E_ref, p0c)
  #r2 = v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]
  #b01 = 2.404825557695773 # first zero of J0
  #d = C_LIGHT*b01/omega
  #arg = (b01*b01)/(d*d)*r2
  #b0, b1 = bessel01_RF(arg)
  #b1 = b1 * b01/d

  t = t0 - v[i,ZI]/C_LIGHT

  #px_0 = v[i,PXI]
  #py_0 = v[i,PYI]
  pz_0 = v[i,PZI]

  phi_particle = omega*t
  #s, c = sincos(phi_particle)

  #coeff = L*E0_over_Rref*b01/(omega*d)*b1*c

  #new_px = px_0 - coeff*v[i,XI]
  #new_py = py_0 - coeff*v[i,YI]
  new_pz = pz_0 + L*E0_over_Rref/C_LIGHT*sin(phi_particle)

  #v[i,PXI] = vifelse(alive, new_px, px_0)
  #v[i,PYI] = vifelse(alive, new_py, py_0)
  v[i,PZI] = vifelse(alive, new_pz, pz_0)

  mad_to_bmad!(i, coords, beta_0, tilde_m, E_ref, p0c)
end


function omega_cavity(i, coords::Coords, a, tilde_m, omega, t0, E0_over_Rref, mm, kn, ks, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)
    #r2 = v[i,XI]*v[i,XI] + v[i,YI]*v[i,YI]
    #b01 = 2.404825557695773 # first zero of J0
    #d = C_LIGHT*b01/omega
    #arg = (b01*b01)/(d*d)*r2
    #b0, b1 = bessel01_RF(arg)
    #b1 = b1 * b01/d
    beta_gamma = (1 + v[i,PZI])/tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    beta = beta_gamma/gamma
    vel = beta*C_LIGHT
    t = t0 - v[i,ZI]/vel

    phi_particle = omega*t
    #s, c = sincos(phi_particle)

    ez = E0_over_Rref*sin(phi_particle)
    ex = zero(ez)
    ey = ex
    e_vec = (ex, ey, ez)

    #coeff = E0_over_Rref/C_LIGHT*b1*c

    bx = ex #-coeff*v[i,YI]
    by = ex #coeff*v[i,XI]
    bz = ex
    b_vec = (bx, by, bz)

    if length(mm) > 0 && mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = ex
      ay = ex
    end

    ox, oy, oz = omega_field(i, coords, a, 0, beta, gamma, ax, ay, e_vec, b_vec, L)
    if length(mm) > 0
      ox1, oy1, oz1 = omega_multipole(i, coords, a, 0, tilde_m, mm, kn, ks, L)
      omega = (ox + ox1, oy + oy1, oz + oz1)
    else
      omega = (ox, oy, oz)
    end
  end end
  return omega
end


"""
This function rotates particle i's quaternion in a cavity.
"""
@makekernel fastgtpsa=true function rotate_spin_cavity!(i, coords::Coords, a, tilde_m, omega, t0, E0_over_Rref, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_cavity(i, coords, a, tilde_m, omega, t0, E0_over_Rref, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end