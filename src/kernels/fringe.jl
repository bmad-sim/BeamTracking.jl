# Curved magnetic
function fringe!(i, coords::Coords, a, tilde_m, beta_0, Kn0, w, w_inv, e1, e2, ::Val{sign}) where {sign}
  @inbounds begin @FastGTPSA begin
    v = coords.v
    rel_p = 1 + v[i,PZI]
    alive = (coords.state[i] == STATE_ALIVE)

    if !isnothing(w)
      rotation!(i, coords, w, 0)
    end

    if sign == 1
      polar_Y!(i, coords, e1, tilde_m, beta_0, a)
    elseif sign == -1
      exact_wedge!(i, coords, -e2, Kn0, tilde_m, beta_0, a)
    end

    if !isnothing(coords.q)
      b_vec = (0, 0, sign*v[i,YI]*Kn0)
      rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
    end
    
    ps2 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
    good_momenta = (ps2 > 0)
    coords.state[i] = vifelse(!good_momenta & alive, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

    xp = v[i,PXI]/ps
    yp = v[i,PYI]/ps
    yp_2 = yp*yp
    yp_factor = 1 + yp_2

    y_sqrt2 = 1 + sign*Kn0*xp*yp/(ps*yp_factor)*2*v[i,YI]
    good_sqrt = (y_sqrt2 > 0)
    coords.state[i] = vifelse(!good_sqrt & alive, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    y_sqrt2_1 = one(y_sqrt2)
    y_sqrt = sqrt(vifelse(good_sqrt, y_sqrt2, y_sqrt2_1))

    new_y  = 2*v[i,YI]/(1 + y_sqrt)
    new_y2 = new_y*new_y
    new_py = v[i,PYI] - sign*Kn0*xp/yp_factor*new_y
    new_x  = v[i,XI]  + sign*Kn0/(2*ps)*new_y2/yp_factor*((1 + xp*xp) - 2*xp*xp*yp_2/yp_factor)
    new_z  = v[i,ZI]  - sign*Kn0*rel_p*xp/(2*ps2)*(1 - yp_2)/(yp_factor*yp_factor)*new_y2

    v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
    v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
    v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
    v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])

    if !isnothing(coords.q)
      b_vec = (0, 0, sign*v[i,YI]*Kn0)
      rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
    end

    if sign == -1
      polar_Y!(i, coords, e2, tilde_m, beta_0, a)
    elseif sign == 1
      exact_wedge!(i, coords, -e1, Kn0, tilde_m, beta_0, a)
    end

    if !isnothing(w_inv)
      rotation!(i, coords, w_inv, 0)
    end
  end end
  return nothing
end


# Straight magnetic
function fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, Ks0, Kn1, Ks1, Kn2, Ks2, ::Val{sign}) where {sign}
  @inbounds begin @FastGTPSA begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)
    rel_p = 1 + v[i,PZI]

    ax = 0
    ay = 0

    if !isnothing(Ksol)
      if sign > 0
        ax = zero(v[i,YI]*Ksol)
        ay = zero(v[i,XI]*Ksol)
      else
        ax = -v[i,YI]*Ksol/2
        ay =  v[i,XI]*Ksol/2
      end
    end

    # Multipole spin
    if !isnothing(coords.q) & (!isnothing(Kn0) | !isnothing(Kn1) | !isnothing(Kn2))
      bz = zero(v[i,XI])
      if !isnothing(Kn0) & !isnothing(Ks0)
        bz += sign*v[i,YI]*Kn0 + sign*v[i,XI]*Ks0
      end
      if !isnothing(Kn1) & !isnothing(Ks1)
        bz += sign*v[i,XI]*v[i,YI]*Kn1 + sign/2*(v[i,XI]*v[i,XI] - v[i,YI]*v[i,YI])*Ks1
      end
      if !isnothing(Kn2) & !isnothing(Ks2)
        x2 = v[i,XI]*v[i,XI]
        y2 = v[i,YI]*v[i,YI]
        bz += sign*v[i,YI]*(x2 - y2/3)*Kn2 + sign*v[i,XI]*(x2/3 - y2)*Ks2
      end
      rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), (0, 0, bz), 1/2)
    end

    # Solenoid spin
    if !isnothing(coords.q) & !isnothing(Ksol)
      b_vec = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, 0)
      rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    end

    # Multipole orbit
    if !isnothing(Kn0) | !isnothing(Kn1) | !isnothing(Kn2)
      fx_over_rel_p = zero(v[i,XI])
      fy_over_rel_p = zero(v[i,XI])
      alphax = zero(v[i,XI])
      alphay = zero(v[i,XI])
      betax = zero(v[i,XI])

      x2 = v[i,XI]*v[i,XI]
      y2 = v[i,YI]*v[i,YI]
      xy = v[i,XI]*v[i,YI]
      x3 = v[i,XI]*x2
      y3 = v[i,YI]*y2
      x4 = x2*x2
      y4 = y2*y2

      if !isnothing(Kn0) & !isnothing(Ks0)
        Kn0_over_rel_p = sign*Kn0/rel_p
        Ks0_over_rel_p = sign*Ks0/rel_p

        fx_over_rel_p += -Kn0_over_rel_p/8*(x2 + 3*y2) - Ks0_over_rel_p/4*xy
        fy_over_rel_p += Kn0_over_rel_p/4*xy + Ks0_over_rel_p/8*(3*x2 + y2)
        alphax += -(Kn0_over_rel_p*v[i,XI] + Ks0_over_rel_p*v[i,YI])/4
        alphay += -(3*Kn0_over_rel_p*v[i,YI] + Ks0_over_rel_p*v[i,XI])/4
        betax += (Kn0_over_rel_p*v[i,YI] + 3*Ks0_over_rel_p*v[i,XI])/4
      end
      if !isnothing(Kn1) & !isnothing(Ks1)
        Kn1_over_rel_p = sign*Kn1/rel_p
        Ks1_over_rel_p = sign*Ks1/rel_p

        fx_over_rel_p += -Kn1_over_rel_p/12*v[i,XI]*(x2 + 3*y2) + Ks1_over_rel_p/6*y3
        fy_over_rel_p += Kn1_over_rel_p/12*v[i,YI]*(3*x2 + y2) + Ks1_over_rel_p/6*x3
        alphax += -Kn1_over_rel_p/4*(x2 + y2)
        alphay += v[i,YI]/2*(Ks1_over_rel_p*v[i,YI] - Kn1_over_rel_p*v[i,XI])
        betax += v[i,XI]/2*(Ks1_over_rel_p*v[i,XI] + Kn1_over_rel_p*v[i,YI])
      end
      if !isnothing(Kn2) & !isnothing(Ks2)
        Kn2_over_rel_p = sign*Kn2/rel_p
        Ks2_over_rel_p = sign*Ks2/rel_p

        fx_over_rel_p += -Kn2_over_rel_p/48*(3*x4 + 6*x2*y2 - 5*y4) + Ks2_over_rel_p/12*xy*(x2 + 3*y2)
        fy_over_rel_p += Kn2_over_rel_p/12*xy*(3*x2 + y2) + Ks2_over_rel_p/48*(5*x4 - 6*x2*y2 - 3*y4)
        alphax += -(Kn2_over_rel_p*v[i,XI] - Ks2_over_rel_p*v[i,YI])/4*(x2 + y2)
        alphay += Kn2_over_rel_p*v[i,YI]/12*(5*y2 - 3*x2) + Ks2_over_rel_p*v[i,XI]/12*(x2 + 9*y2)
        betax += Kn2_over_rel_p*v[i,YI]/12*(y2 + 9*x2) + Ks2_over_rel_p*v[i,XI]/12*(5*x2 - 3*y2)
      end

      delta = 1 - alphax*alphax - alphay*betax
      good_delta = (delta > 0)
      alive_at_start = (coords.state[i] == STATE_ALIVE)
      coords.state[i] = vifelse(!good_delta & alive_at_start, STATE_LOST, coords.state[i])
      alive = (coords.state[i] == STATE_ALIVE)
      delta_1 = one(delta)
      delta = vifelse(good_delta, delta, delta_1)
      px_over_delta = v[i,PXI]/delta
      py_over_delta = v[i,PYI]/delta

      new_x  = v[i,XI] - fx_over_rel_p
      new_y  = v[i,YI] - fy_over_rel_p
      new_px = (1 + alphax)*px_over_delta + betax*py_over_delta
      new_py = alphay*px_over_delta + (1 - alphax)*py_over_delta
      new_z  = v[i,ZI] + new_px/rel_p*fx_over_rel_p + new_py/rel_p*fy_over_rel_p

      v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
      v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
      v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
      v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
      v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])
  end

    # Solenoid orbit
    if !isnothing(Ksol)
      if sign > 0
        ax = -v[i,YI]*Ksol/2
        ay =  v[i,XI]*Ksol/2
      else
        ax = zero(v[i,YI]*Ksol)
        ay = zero(v[i,XI]*Ksol)
      end
    end

    # Multipole spin
    if !isnothing(coords.q)
      bz = zero(v[i,XI])
      if !isnothing(Kn0) & !isnothing(Ks0)
        bz += sign*v[i,YI]*Kn0 + sign*v[i,XI]*Ks0
      end
      if !isnothing(Kn1) & !isnothing(Ks1)
        bz += sign*v[i,XI]*v[i,YI]*Kn1 + sign/2*(v[i,XI]*v[i,XI] - v[i,YI]*v[i,YI])*Ks1
      end
      if !isnothing(Kn2) & !isnothing(Ks2)
        x2 = v[i,XI]*v[i,XI]
        y2 = v[i,YI]*v[i,YI]
        bz += sign*v[i,YI]*(x2 - y2/3)*Kn2 + sign*v[i,XI]*(x2/3 - y2)*Ks2
      end
      rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), (0, 0, bz), 1/2)
    end

    # Solenoid spin
    if !isnothing(coords.q) & !isnothing(Ksol)
      b_vec = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, 0)
      rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    end
  end end
  return nothing
end


# Straight electric
function fringe!(i, coords::Coords, a, tilde_m, kE, w, w_inv, ::Val{sign}) where {sign}
  @inbounds begin @FastGTPSA begin
    if !isnothing(coords.q)
      if !isnothing(w)
        rotation!(i, coords, w, 0)
      end

      beta_0 = 1/sqrt(1 + tilde_m*tilde_m)
      phi = -kE*coords.v[i,XI]
      e_vec = (0, 0, -sign*c_light(eltype(coords.v))*phi)

      if sign > 0
        phi_in = zero(phi)
        phi_out = phi
      else
        phi_in = phi
        phi_out = zero(phi)
      end

      mad_to_bmad!(i, coords, beta_0, tilde_m, phi_in)
      rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, e_vec, (0, 0, 0), 1/2)
      bmad_to_mad!(i, coords, beta_0, tilde_m, phi_in)
      mad_to_bmad!(i, coords, beta_0, tilde_m, phi_out)
      rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, e_vec, (0, 0, 0), 1/2)
      bmad_to_mad!(i, coords, beta_0, tilde_m, phi_out)

      if !isnothing(w_inv)
        rotation!(i, coords, w_inv, 0)
      end
    end
  end end
  return nothing
end


@makekernel fastgtpsa=true function exact_wedge!(i, coords::Coords, theta, Kn0, tilde_m, beta_0, a)
  v = coords.v
  px_0 = v[i,PXI]
  py_0 = v[i,PYI]
  rel_p = 1 + v[i,PZI]
  rel_p2 = rel_p*rel_p

  pt2 = rel_p2 - v[i,PYI]*v[i,PYI]
  good_momenta = (pt2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  pt2_1 = one(pt2)
  pt = sqrt(vifelse(good_momenta, pt2, pt2_1))

  arg = v[i,PXI]/pt
  abs_arg = abs(arg)
  arg_1 = one(arg)
  arg_0 = zero(arg)
  good_arg = (abs_arg < arg_1)
  coords.state[i] = vifelse(!good_arg & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)

  phi1 = theta + asin(vifelse(good_arg, arg, arg_0))
  gp = Kn0/pt
  splus, cplus = sincos(phi1)
  sth, cth = sincos(theta)
  sgn = vifelse(theta > 0, one(theta), -one(theta))
  s = v[i,XI]*sth

  cond = cplus*cplus + gp*s*(2*splus - gp*s)
  good_cond = (cond > 0)
  coords.state[i] = vifelse(!good_cond & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  cond_1 = one(cond)
  nasty_sqrt = sqrt(vifelse(good_cond, cond, cond_1))

  gp_0 = zero(gp)
  abs_gp = abs(gp)
  good_gp = (abs_gp > gp_0)
  gp_1 = one(gp)
  gp_safe = vifelse(good_gp, gp, gp_1)
  pos_cplus = (cplus > 0)
  branch1 = !good_gp | pos_cplus
  den = vifelse(branch1, nasty_sqrt + cplus, one(cplus))
  beta = (2*splus - gp*s)/den
  xi1 = s*beta
  xi2 = (nasty_sqrt - cplus)/gp_safe
  xi = vifelse(branch1, xi1, xi2)

  u = sgn*s
  Lcv = -u
  e = vifelse(u > 0, one(u), -one(u))  
  ang1 = sgn*atan2(sgn*e*beta, e)          
  ang2 = sgn*atan2(xi, u)                   
  ang = vifelse(branch1, ang1, ang2)
  thetap = 2*(phi1 - ang)
  Lp = sgn*sqrt(Lcv*Lcv + xi*xi)/sincu(thetap/2) 

  new_x = v[i,XI]*cth + xi
  new_px = pt*sin(phi1 - thetap)
  new_y = v[i,YI] + v[i,PYI]*Lp/pt
  new_z = v[i,ZI] - rel_p*Lp/pt
  v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    ps2 = pt2 - px_0*px_0
    good_momenta = (ps2 > 0)
    alive_at_start = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

    beta_gamma = rel_p/tilde_m
    beta_gamma2 = beta_gamma*beta_gamma
    gamma_minus_1 = beta_gamma2/(1 + sqrt(1 + beta_gamma2))
    gamma = gamma_minus_1 + 1
    vt_over_rel_p = Lp/pt
    factor = vt_over_rel_p*Kn0
    coeff = a*factor*gamma_minus_1*py_0/rel_p2

    o1 = coeff*px_0
    o2 = coeff*py_0 - a*gamma*factor
    o3 = coeff*ps
    q1 = expq((o1, o2, o3), alive)

    theta_0 = zero(theta - factor)
    angle = vifelse(alive, (theta - factor)/2, theta_0)
    s, c = sincos(angle)
    q2 = (c, theta_0, s, theta_0)
    q3 = quat_mul(q2, q1)
    q4 = quat_mul(q3, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q4
  end
end


@makekernel fastgtpsa=true function polar_Y!(i, coords::Coords, theta, tilde_m, beta_0, a)
  v = coords.v
  rel_p = 1 + v[i,PZI]

  ps2 = rel_p*rel_p - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI]
  good_momenta = (ps2 > 0)
  alive_at_start = (coords.state[i] == STATE_ALIVE)
  coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  ps2_1 = one(ps2)
  ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

  s, c = sincos(theta)
  t = s/c

  denom = 1 - t*v[i,PXI]/ps
  abs_denom = abs(denom)
  denom_1 = one(denom)
  denom_0 = zero(denom)
  good_denom = (abs_denom > denom_0)
  coords.state[i] = vifelse(!good_denom & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  denom = vifelse(good_denom, denom, denom_1)

  new_x = v[i,XI]/(c*denom) 
  coeff = s*new_x/ps
  new_px = v[i,PXI]*c + ps*s
  new_y = v[i,YI] + v[i,PYI]*coeff
  new_z = v[i,ZI] - rel_p*coeff
  v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
  v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])

  if !isnothing(coords.q) && !isnothing(a)
    q = coords.q

    theta_0 = zero(theta)
    angle = vifelse(alive, theta/2, theta_0)
    s, c = sincos(angle)
    q1 = (c, theta_0, s, theta_0)
    q2 = quat_mul(q1, q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
    q[i,Q0], q[i,QX], q[i,QY], q[i,QZ] = q2
  end
end