# Curved magnetic
#
# Second order hard edge bend fringe of Hwang and Lee with the quadrupole terms of
# Iselin (Bmad manual "Bend Second Order Fringe Map"). The Lie generator is
#   Ω = K(x, y, pz) + sign * B(x, px, y, py, pz)
# with K independent of the transverse momenta and B linear in them. Instead of
# truncating exp(:Ω:) at second order, which is not symplectic, the map is
#   exp(:K/2:) exp(:sign*B1:) exp(:sign*B2:) exp(:sign*B3:) exp(:K/2:)
# where B = B1 + B2 + B3 and each factor is integrated exactly. Since {K, {K, B}} = 0,
# this agrees with exp(:Ω:) through second order and is exactly symplectic.
# edge1_int and edge2_int are fint*hgap in Bmad.
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Kn0, Kn1, w, w_inv, e1, e2, edge1_int, edge2_int, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  if sign > 0
    e = e1
    edge_int = edge1_int
  else
    e = e2
    edge_int = edge2_int
  end

  if isnothing(Kn1)
    k1 = zero(Kn0)
  else
    k1 = Kn1
  end

  sn, cs = sincos(e)
  t = sn/cs
  t2 = t*t
  sec2 = 1/(cs*cs)
  f = Kn0*t

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  if !isnothing(coords.q)
    b_vec = (-v[i,YI]*f, -v[i,XI]*f, sign*v[i,YI]*Kn0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
  end

  rel_p = 1 + v[i,PZI]
  # K = f*(x^2 - y^2)/2 + (cy2*y^2 + cx3*x^3 + cxy2*x*y^2)/rel_p
  cy2  = Kn0*Kn0*sec2/cs*(1 + sn*sn)*edge_int
  cx3  = (4*k1*t - Kn0*Kn0*t2*t)/12
  cxy2 = (-4*k1*t + Kn0*Kn0*t*sec2)/4
  # B = c*(t2*x^2*px - 2*t2*x*y*py - sec2*y^2*px)
  c = sign*Kn0/(2*rel_p)

  x  = v[i,XI]
  px = v[i,PXI]
  y  = v[i,YI]
  py = v[i,PYI]
  z  = v[i,ZI]

  # exp(:K/2:)
  px = px + (f*x + (3*cx3*x*x + cxy2*y*y)/rel_p)/2
  py = py + (-f*y + 2*(cy2 + cxy2*x)*y/rel_p)/2
  z  = z  + (cy2*y*y + cx3*x*x*x + cxy2*x*y*y)/(2*rel_p*rel_p)

  # exp(:B1:), B1 = c*t2*x^2*px
  u = 1 + c*t2*x
  good_u = (u > 0)
  coords.state[i] = vifelse(!good_u & alive, STATE_LOST, coords.state[i])
  alive = (coords.state[i] == STATE_ALIVE)
  u = vifelse(good_u, u, one(u))
  z  = z + c*t2*x*x*px/rel_p
  x  = x/u
  px = px*u*u

  # exp(:B2:), B2 = -c*sec2*y^2*px
  z  = z - c*sec2*y*y*px/rel_p
  x  = x + c*sec2*y*y
  py = py - 2*c*sec2*y*px

  # exp(:B3:), B3 = -2*c*t2*x*y*py
  z  = z - 2*c*t2*x*y*py/rel_p
  px = px - 2*c*t2*y*py
  ex = exp(2*c*t2*x)
  y  = y*ex
  py = py/ex

  # exp(:K/2:)
  px = px + (f*x + (3*cx3*x*x + cxy2*y*y)/rel_p)/2
  py = py + (-f*y + 2*(cy2 + cxy2*x)*y/rel_p)/2
  z  = z  + (cy2*y*y + cx3*x*x*x + cxy2*x*y*y)/(2*rel_p*rel_p)

  v[i,XI]  = vifelse(alive, x,  v[i,XI])
  v[i,PXI] = vifelse(alive, px, v[i,PXI])
  v[i,YI]  = vifelse(alive, y,  v[i,YI])
  v[i,PYI] = vifelse(alive, py, v[i,PYI])
  v[i,ZI]  = vifelse(alive, z,  v[i,ZI])

  if !isnothing(coords.q)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
  end

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end
end


# Curved magnetic, linear hard edge kick only. Used by exact bend tracking.
@makekernel fastgtpsa=true function linear_bend_fringe!(i, coords::Coords, a, tilde_m, Kn0, w, w_inv, e1, e2, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  if sign > 0
    e = e1
  else
    e = e2
  end

  f = Kn0*tan(e)

  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  if !isnothing(coords.q)
    b_vec = (-v[i,YI]*f, -v[i,XI]*f, sign*v[i,YI]*Kn0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
  end

  new_px = v[i,PXI] + f*v[i,XI]
  new_py = v[i,PYI] - f*v[i,YI]
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])

  if !isnothing(coords.q)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
  end

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end
end


# Straight magnetic
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, w0, w0_inv, Kn1, w1, w1_inv, sign)
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

  # Quadrupole
  if !isnothing(Kn1) && !isnothing(coords.q)
    if !isnothing(w1)
      rotation!(i, coords, w1, 0)
    end
    b_vec = (0, 0, sign*v[i,XI]*v[i,YI]*Kn1)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    if !isnothing(w1_inv)
      rotation!(i, coords, w1_inv, 0)
    end
  end

  # Dipole
  if !isnothing(Kn0) && !isnothing(coords.q)
    if !isnothing(w0)
      rotation!(i, coords, w0, 0)
    end
    b_vec = (0, 0, sign*v[i,YI]*Kn0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    if !isnothing(w0_inv)
      rotation!(i, coords, w0_inv, 0)
    end
  end

  # Solenoid
  if !isnothing(Ksol) && !isnothing(coords.q)
    b_vec = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, 0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
  end

  # Quadrupole
  if !isnothing(Kn1)
    if !isnothing(w1)
      rotation!(i, coords, w1, 0)
    end

    Kn1_over_rel_p = Kn1/rel_p

    x2 = v[i,XI]*v[i,XI]
    x3 = v[i,XI]*x2
    y2 = v[i,YI]*v[i,YI]
    y3 = v[i,YI]*y2
    x2y = x2*v[i,YI]
    y2x = y2*v[i,XI]

    alphax = -sign*Kn1_over_rel_p/4*(x2 + y2)
    alphay = -sign*Kn1_over_rel_p/2*v[i,XI]*v[i,YI]
    delta = 1 - alphax*alphax + alphay*alphay
    px_over_delta = v[i,PXI]/delta
    py_over_delta = v[i,PYI]/delta

    new_x  = v[i,XI] + sign*Kn1_over_rel_p/12*(x3 + 3*y2x)
    new_y  = v[i,YI] - sign*Kn1_over_rel_p/12*(y3 + 3*x2y)
    new_px = (1 + alphax)*px_over_delta -       alphay*py_over_delta
    new_py =       alphay*px_over_delta + (1 - alphax)*py_over_delta
    new_z  = v[i,ZI] + sign*Kn1_over_rel_p/(12*rel_p)*(y3*new_py - x3*new_px + 3*x2y*new_py - 3*y2x*new_px)

    v[i,XI]  = vifelse(alive, new_x,  v[i,XI])
    v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
    v[i,YI]  = vifelse(alive, new_y,  v[i,YI])
    v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
    v[i,ZI]  = vifelse(alive, new_z,  v[i,ZI])

    if !isnothing(w1_inv)
      rotation!(i, coords, w1_inv, 0)
    end
  end

  # Dipole
  if !isnothing(Kn0)
    if !isnothing(w0)
      rotation!(i, coords, w0, 0)
    end

    px = v[i,PXI] - ax
    py = v[i,PYI] - ay
    ps2 = rel_p*rel_p - px*px - py*py
    good_momenta = (ps2 > 0)
    coords.state[i] = vifelse(!good_momenta & alive, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

    xp = px/ps
    yp = py/ps
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

    if !isnothing(w0_inv)
      rotation!(i, coords, w0_inv, 0)
    end
  end

  # Solenoid
  if !isnothing(Ksol)
    if sign > 0
      ax = -v[i,YI]*Ksol/2
      ay =  v[i,XI]*Ksol/2
    else
      ax = zero(v[i,YI]*Ksol)
      ay = zero(v[i,XI]*Ksol)
    end
  end

  # Quadrupole
  if !isnothing(Kn1) && !isnothing(coords.q)
    if !isnothing(w1)
      rotation!(i, coords, w1, 0)
    end
    b_vec = (0, 0, sign*v[i,XI]*v[i,YI]*Kn1)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    if !isnothing(w1_inv)
      rotation!(i, coords, w1_inv, 0)
    end
  end

  # Dipole
  if !isnothing(Kn0) && !isnothing(coords.q)
    if !isnothing(w0)
      rotation!(i, coords, w0, 0)
    end
    b_vec = (0, 0, sign*v[i,YI]*Kn0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
    if !isnothing(w0_inv)
      rotation!(i, coords, w0_inv, 0)
    end
  end

  # Solenoid
  if !isnothing(Ksol) && !isnothing(coords.q)
    b_vec = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, 0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
  end
end


# Straight electric
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, kE, w, w_inv, sign)
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
end