# Curved
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Kn0, w, w_inv, e, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

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
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, Ks0, Kn1, w, w_inv, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)

  rel_p = 1 + v[i,PZI]
  Kn1_over_rel_p = Kn1/rel_p

  if !isnothing(coords.q)
    s = 2*w[1]*w[4]
    c = w[1]^2 - w[4]^2
    b = (-sign*v[i,XI]*Ksol/2, -sign*v[i,YI]*Ksol/2, sign*v[i,YI]*Kn0 + sign*v[i,XI]*Ks0)
    b_rotated = (b[1]*c + b[2]*s, b[2]*c - b[1]*s, b[3]) # rotate opposite sense of w
  end

  # Need to rotate so that quadrupole component is upright
  if !isnothing(w)
    rotation!(i, coords, w, 0)
  end

  if !isnothing(coords.q)
    b_quad = (0, 0, 0)
    b_vec = b_rotated .+ b_quad
    rotate_spin_field!(i, coords, a, 0, tilde_m, 0, 0, (0, 0, 0), b_vec, 1/2)
  end

  x2 = v[i,XI]*v[i,XI]
  x3 = x2*v[i,XI]
  y2 = v[i,YI]*v[i,YI]
  y3 = y2*v[i,YI]
  new_x = v[i,XI] + sign*Kn1_over_rel_p/12*(x3 + 3*y2*v[i,XI])
  new_px = v[i,PXI] + sign*Kn1_over_rel_p/4*(2*v[i,XI]*v[i,YI]*v[i,PYI] - (x2 + y2)*v[i,PXI]) 
  new_y = v[i,YI] - sign*Kn1_over_rel_p/12*(y3 + 3*x2*v[i,YI])
  new_py = v[i,PYI] - sign*Kn1_over_rel_p/4*(2*v[i,XI]*v[i,YI]*v[i,PXI] - (x2 + y2)*v[i,PYI])
  new_z = v[i,ZI] + sign*Kn1_over_rel_p/12*(y3*v[i,PYI] - x3*v[i,PXI] + 3*x2*v[i,YI]*v[i,PYI] - 3*y2*v[i,XI]*v[i,PXI])
  v[i,XI] =  vifelse(alive, new_x,  v[i,XI])
  v[i,PXI] = vifelse(alive, new_px, v[i,PXI])
  v[i,YI] =  vifelse(alive, new_y,  v[i,YI])
  v[i,PYI] = vifelse(alive, new_py, v[i,PYI])
  v[i,ZI] =  vifelse(alive, new_z,  v[i,ZI])

  if !isnothing(coords.q)
    ax = -v[i,YI]*Ksol/2
    ay =  v[i,XI]*Ksol/2
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
  end

  if !isnothing(w_inv)
    rotation!(i, coords, w_inv, 0)
  end
end