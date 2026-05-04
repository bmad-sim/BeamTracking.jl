@inline function cavity!(i, coords::Coords, s, radiation_params, beta_0, gamsqr_0, tilde_m, E_ref, a, omega, t_ref, E0_normalized, Ksol, ::Val{sol}, mm, kn, ks, L) where {sol}
  @inbounds begin @FastGTPSA begin
    #t_ref += (s + L/2)/(beta_0*C_LIGHT)
    multipoles = (length(mm) > 0)
    if sol
      exact_solenoid!(i, coords, s, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
    else
      exact_drift!(i, coords, s, beta_0, gamsqr_0, tilde_m, L / 2)
    end

    if !isnothing(radiation_params)
      q, mc2, E_ref = radiation_params
      deterministic_radiation_cavity!(i, coords, q, mc2, E_ref, omega, t_ref, E0_normalized, mm, kn, ks, L / 2)
    end

    if multipoles
      knl = kn .* L ./ 2
      ksl = ks .* L ./ 2
      multipole_kick!(i, coords, mm, knl, ksl, -1)
    end

    if isnothing(coords.q)
      cavity_kick!(i, coords, beta_0, tilde_m, omega, t_ref, E0_normalized, L)
    else
      cavity_kick!(i, coords, beta_0, tilde_m, omega, t_ref, E0_normalized, L / 2)
      rotate_spin_cavity!(i, coords, a, tilde_m, omega, t_ref, E0_normalized, mm, kn, ks, L)
      cavity_kick!(i, coords, beta_0, tilde_m, omega, t_ref, E0_normalized, L / 2)
    end

    if multipoles
      multipole_kick!(i, coords, mm, knl, ksl, -1)
    end

    if !isnothing(radiation_params)
      deterministic_radiation_cavity!(i, coords, q, mc2, E_ref, omega, t_ref, E0_normalized, mm, kn, ks, L / 2)
    end

    if sol
      exact_solenoid!(i, coords, s, Ksol, beta_0, gamsqr_0, tilde_m, L / 2)
    else
      exact_drift!(i, coords, s, beta_0, gamsqr_0, tilde_m, L / 2)
    end
  end end
  return nothing
end


"""
Converts the longitudinal coordinates from (z, pz) to (τ, pτ) where
τ = c(t_ref - t) and pτ = (E - E_ref)/pc_ref. Making this well-conditioned is
nontrivial and this implementation still may not be optimal.
"""
@makekernel fastgtpsa=true function bmad_to_mad!(i, coords::Coords, beta_0, tilde_m, phi)
  v = coords.v

  pz = v[i,PZI]
  rel_p = 1 + pz
  y = beta_0*(2*pz + pz*pz)
  good = (beta_0*y > -1)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  y_1 = one(y)
  x = 1 + beta_0*y
  ptau = y/(1 + sqrt(vifelse(good, x, y_1))) + phi
  beta = rel_p/sqrt(rel_p*rel_p + tilde_m*tilde_m)
  tau = v[i,ZI]/beta

  v[i,ZI]  = vifelse(alive, tau, v[i,ZI])
  v[i,PZI] = vifelse(alive, ptau, v[i,PZI])
end


"""
Converts the longitudinal coordinates from (τ, pτ) to (z, pz) where
τ = c(t_ref - t) and pτ = (E - E_ref)/pc_ref. Making this well-conditioned is
nontrivial and this implementation still may not be optimal.
"""
@makekernel fastgtpsa=true function mad_to_bmad!(i, coords::Coords, beta_0, tilde_m, phi)
  v = coords.v

  ptau = v[i,PZI]
  y = ptau*2/beta_0 + ptau*ptau - phi*2/beta_0 + phi*phi - 2*ptau*phi
  good = (y > -1)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  y_1 = one(y)
  x = 1 + y
  pz = y/(1 + sqrt(vifelse(good, x, y_1)))
  rel_p = 1 + pz
  beta = rel_p/sqrt(rel_p*rel_p + tilde_m*tilde_m)
  z = v[i,ZI]*beta

  v[i,ZI]  = vifelse(alive, z, v[i,ZI])
  v[i,PZI] = vifelse(alive, pz, v[i,PZI])
end


@makekernel fastgtpsa=true function cavity_kick!(i, coords::Coords, beta_0, tilde_m, omega, t_ref, E0_normalized, L)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  bmad_to_mad!(i, coords, beta_0, tilde_m, 0)

  t = t_ref - v[i,ZI]/C_LIGHT
  new_pz = v[i,PZI] + L*E0_normalized/C_LIGHT*sin(omega*t)

  v[i,PZI] = vifelse(alive, new_pz, v[i,PZI])
  mad_to_bmad!(i, coords, beta_0, tilde_m, 0)
end


"""
Returns the integrated spin-precession vector for an RF cavity, possibly with
multipoles.
"""
function omega_cavity(i, coords::Coords, a, tilde_m, omega, t_ref, E0_normalized, mm, kn, ks, L)
  @FastGTPSA begin @inbounds begin
    v = coords.v

    beta_gamma = (1 + v[i,PZI])/tilde_m
    gamma = sqrt(1 + beta_gamma*beta_gamma)
    beta = beta_gamma/gamma
    vel = beta*C_LIGHT
    t = t_ref - v[i,ZI]/vel

    ez = E0_normalized*sin(omega*t)
    ex = zero(ez)
    ey = ex
    e_vec = (ex, ey, ez)
    b_vec = (ex, ex, ex)
    if length(mm) > 0 && mm[1] == 0
      ax = -v[i,YI] * kn[1] / 2
      ay =  v[i,XI] * kn[1] / 2
    else
      ax = ex
      ay = ex
    end

    ox, oy, oz = omega_field(i, coords, a, 0, tilde_m, ax, ay, e_vec, b_vec, L)
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
Gives radiation damping kick in an RF cavity, possibly with multipoles.
"""
@makekernel fastgtpsa=true function deterministic_radiation_cavity!(i, coords::Coords, q, mc2, E_ref, omega, t_ref, E0_normalized, mm, kn, ks, L) 
  v = coords.v

  t = t_ref - v[i,ZI]/C_LIGHT # ultrarelativistic radiation
  ez = E0_normalized*sin(omega*t)
  ex = zero(ez)
  ey = ex
  e_vec = (ex, ey, ez)
  
  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
    b_vec = (bx, by, kn[1])
  else
    ax = ex
    ay = ex
    b_vec = (bx, by, ex)
  end

  deterministic_radiation_field!(i, coords, q, mc2, E_ref, 0, ax, ay, e_vec, b_vec, L)
end



"""
Gives radiation diffusion kick in an RF cavity, possibly with multipoles.
"""
@makekernel fastgtpsa=true function stochastic_radiation!(i, coords::Coords, s, ::typeof(cavity!), backend, q, mc2, E_ref, omega, t_ref, E0_normalized, mm, kn, ks, L) 
  v = coords.v

  t = t_ref - v[i,ZI]/C_LIGHT # ultrarelativistic radiation
  ez = E0_normalized*sin(omega*t)
  ex = zero(ez)
  ey = ex
  e_vec = (ex, ey, ez)
  
  bx, by = normalized_field(mm, kn, ks, v[i,XI], v[i,YI], -1)
  if mm[1] == 0
    ax = -v[i,YI] * kn[1] / 2
    ay =  v[i,XI] * kn[1] / 2
    b_vec = (bx, by, kn[1])
  else
    ax = ex
    ay = ex
    b_vec = (bx, by, ex)
  end

  stochastic_radiation_field!(i, coords, backend, q, mc2, E_ref, 0, ax, ay, e_vec, b_vec, L)
end


"""
This function rotates particle i's quaternion in a cavity.
"""
@makekernel fastgtpsa=true function rotate_spin_cavity!(i, coords::Coords, a, tilde_m, omega, t_ref, E0_normalized, mm, kn, ks, L)
  q2 = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  q1 = expq(omega_cavity(i, coords, a, tilde_m, omega, t_ref, E0_normalized, mm, kn, ks, L), alive)
  q3 = quat_mul(q1, q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ])
  q2[i,Q0], q2[i,QX], q2[i,QY], q2[i,QZ] = q3
end