# Curved
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Kn0, w, w_inv, e1, e2, sign)
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


# Straight
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, w0, w0_inv, Kn1, w1, w1_inv, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p = 1 + v[i,PZI]

  # Quadrupole
  if !isnothing(Kn1) && !isnothing(coords.q)
    if !isnothing(w1)
      rotation!(i, coords, w1, 0)
    end
    b_vec = (0, 0, sign*v[i,XI]*v[i,YI]*Kn1)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
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
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
    if !isnothing(w0_inv)
      rotation!(i, coords, w0_inv, 0)
    end
  end

  # Solenoid
  if !isnothing(Ksol) && !isnothing(coords.q)
    b_vec = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, 0)
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
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

    alphax = -sign*Kn1_over_rel_p/4*(x2 + y2)
    alphay = -sign*Kn1_over_rel_p/2*v[i,XI]*v[i,YI]
    delta = 1 - alphax*alphax + alphay*alphay
    px_over_delta = v[i,PXI]/delta
    py_over_delta = v[i,PYI]/delta

    new_x  = v[i,XI] + sign*Kn1_over_rel_p/12*(x3 + 3*y2*v[i,XI])
    new_y  = v[i,YI] - sign*Kn1_over_rel_p/12*(y3 + 3*x2*v[i,YI])
    new_px = (1 + alphax)*px_over_delta -       alphay*py_over_delta
    new_py =       alphay*px_over_delta + (1 - alphax)*py_over_delta
    new_z  = v[i,ZI] + sign*Kn1_over_rel_p/(12*rel_p)*(y3*new_py - x3*new_px + 3*x2*v[i,YI]*new_py - 3*y2*v[i,XI]*new_px)

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

    px = v[i,PXI]
    py = v[i,PYI]
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
    ax = -v[i,YI]*Ksol/2
    ay =  v[i,XI]*Ksol/2
  else
    ax = 0
    ay = 0
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