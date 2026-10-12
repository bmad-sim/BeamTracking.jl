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
@makekernel fastgtpsa=true function fringe!(i, coords::Coords, a, tilde_m, Ksol, Kn0, w0, w0_inv, mm, kn, ks, sign)
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

  # Quadrupole and higher multipoles
  if !isnothing(mm) && !isnothing(coords.q)
    b_vec = (0, 0, sign*multipole_fringe_bs(mm, kn, ks, v[i,XI], v[i,YI]))
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
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

  # Quadrupole and higher multipoles
  if !isnothing(mm)
    multipole_fringe!(i, coords, mm, kn, ks, sign)
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

  # Quadrupole and higher multipoles
  if !isnothing(mm) && !isnothing(coords.q)
    b_vec = (0, 0, sign*multipole_fringe_bs(mm, kn, ks, v[i,XI], v[i,YI]))
    rotate_spin_field!(i, coords, a, 0, tilde_m, ax, ay, (0, 0, 0), b_vec, 1/2)
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

# =========== Multipole hard edge fringe (quadrupole and higher) ============= #
#
# Hard edge fringe of J. S. Berg, "Higher Order Hard Edge End Field Effects", EPAC 2004,
# for straight multipoles of order m >= 2 (m = 2 is a quadrupole) with the Maxwell
# consistent field
#   ψ = Σ_j c_j S^(2j)(s) r^(2j) Im(K w^m)/m!,  c_j = (-1)^j m!/(4^j j! (j+m)!),
# where w = x + i y, K = Kn + i Ks, and S steps from 0 to 1. In the gauge with
# A_s = -S Re(K w^m)/m! (the body gauge) and A_⊥ the rotated gradient of
#   χ = -Σ_j c_(j+1) S^(2j+1)(s) r^(2j+2) Im(K w^m)/m!,
# only the transverse vector potential contributes and the entrance generator is
#   W = Σ_(j=0..2) Re(P̄ D^(2j) G_j),  G_j = 2i ∂χ_j/∂w̄ (without S),
# with P = x' + i y' = (px + i py)/ps and D = P ∂/∂w + P̄ ∂/∂w̄. This keeps terms through
# fifth order in the angles. Collecting terms, W = Re(P̄ (E A + Ē B)) with E = K w^(m-1),
# Ē = conj(K) w̄^(m-2) and A, B low order polynomials in w, w̄, P, P̄ given below.
# The exit generator is -W. The map z_f = z_i + J ∇W((z_i + z_f)/2) is solved with the
# implicit midpoint rule, which is symplectic. W is linear in Kn and Ks so the map is
# differentiable everywhere, including at zero strength.

const MULTIPOLE_FRINGE_ITERATIONS = 6

# Complex arithmetic on (re, im) tuples, usable with any number type (TPS, SIMD.Vec, Dual).
@inline cmul(a, b) = (a[1]*b[1] - a[2]*b[2], a[1]*b[2] + a[2]*b[1])
@inline cadd(a, b) = (a[1] + b[1], a[2] + b[2])
@inline cscale(c, a) = (c*a[1], c*a[2])
@inline cconj(a) = (a[1], -a[2])
@inline function cpow(a, n)
  r = (one(a[1]), zero(a[1]))
  for _ in 1:n
    r = cmul(r, a)
  end
  return r
end

"""
    multipole_fringe_gradient(m, kn, ks, x, y, xp, yp)

Gradient (∂W/∂x, ∂W/∂y, ∂W/∂x', ∂W/∂y') of the entrance hard edge fringe generator W
of a multipole of order `m >= 2` with normal and skew strengths `kn`, `ks`, where
x' = px/ps and y' = py/ps.
"""
@inline function multipole_fringe_gradient(m, kn, ks, x, y, xp, yp)
  fac1 = 1.0
  for k in 2:m+1
    fac1 *= k
  end
  fac2 = fac1*(m+2)
  fac3 = fac2*(m+3)
  p0 = -1/(4*fac1)
  p1 =  1/(32*fac2)
  p2 = -1/(384*fac3)

  q1a = 2*p1*(m+2)*(m+1)
  q1b = 4*p1*(m+2)
  q2a = 36*p2*(m+3)*(m+2)
  q2b = 24*p2*(m+3)*(m+2)*(m+1)
  q2c = 3*p2*(m+3)*(m+2)*(m+1)*m
  r0  = -p0*(m+1)
  r1a = -2*p1*(m+2)
  r1b = -4*p1*(m+2)*(m+1)
  r1c = -p1*(m+2)*m*(m+1)
  r2a = -p2*(m+3)*(m+2)*(m+1)*m*(m-1)
  r2b = -12*p2*(m+3)*(m+2)*(m+1)*m
  r2c = -36*p2*(m+3)*(m+2)*(m+1)
  r2d = -24*p2*(m+3)*(m+2)

  w   = (x, y)
  wb  = (x, -y)
  P   = (xp, yp)
  Pb  = (xp, -yp)
  r2  = x*x + y*y
  pp2 = xp*xp + yp*yp
  w2  = cmul(w, w)
  w3  = cmul(w2, w)
  wb2 = cconj(w2)
  wb3 = cconj(w3)
  P2  = cmul(P, P)
  P3  = cmul(P2, P)
  P4  = cmul(P2, P2)
  Pb2 = cconj(P2)
  Pb3 = cconj(P3)
  Pb4 = cconj(P4)

  # A = a1 w² + a2 r² P² + a3 w̄² P⁴,  B = b1 r² w̄ + b2 w̄³ P² + b3 r² w P̄² + b4 w³ P̄⁴
  # with a_i, b_i functions of |P|².
  a1 = p0 + q1b*pp2 + q2a*pp2*pp2
  a2 = q1a + q2b*pp2
  a3 = q2c
  b1 = r0 + r1b*pp2 + r2c*pp2*pp2
  b2 = r1a + r2d*pp2
  b3 = r1c + r2b*pp2
  b4 = r2a
  da1 = q1b + 2*q2a*pp2
  da2 = q2b
  db1 = r1b + 2*r2c*pp2
  db2 = r2d
  db3 = r2b

  r2P2   = cscale(r2, P2)
  wbP4   = cmul(wb, P4)
  wb2P4  = cmul(wb2, P4)
  r2wb   = cscale(r2, wb)
  wb3P2  = cmul(wb3, P2)
  r2wPb2 = cscale(r2, cmul(w, Pb2))
  w3Pb4  = cmul(w3, Pb4)

  A    = cadd(cadd(cscale(a1, w2), cscale(a2, r2P2)), cscale(a3, wb2P4))
  A_w  = cadd(cscale(2*a1, w), cscale(a2, cmul(wb, P2)))
  A_wb = cadd(cscale(a2, cmul(w, P2)), cscale(2*a3, wbP4))
  A_P  = cadd(cadd(cmul(Pb, cadd(cscale(da1, w2), cscale(da2, r2P2))), cscale(2*a2*r2, P)), cscale(4*a3, cmul(wb2, P3)))
  A_Pb = cmul(P, cadd(cscale(da1, w2), cscale(da2, r2P2)))

  B    = cadd(cadd(cscale(b1, r2wb), cscale(b2, wb3P2)), cadd(cscale(b3, r2wPb2), cscale(b4, w3Pb4)))
  B_w  = cadd(cadd(cscale(b1, wb2), cscale(2*b3*r2, Pb2)), cscale(3*b4, cmul(w2, Pb4)))
  B_wb = cadd(cadd(cscale(2*b1*r2, (one(x), zero(x))), cscale(3*b2, cmul(wb2, P2))), cscale(b3, cmul(w2, Pb2)))
  B_P  = cadd(cmul(Pb, cadd(cadd(cscale(db1, r2wb), cscale(db2, wb3P2)), cscale(db3, r2wPb2))), cscale(2*b2, cmul(wb3, P)))
  B_Pb = cadd(cadd(cmul(P, cadd(cadd(cscale(db1, r2wb), cscale(db2, wb3P2)), cscale(db3, r2wPb2))),
                   cscale(2*b3*r2, cmul(w, Pb))), cscale(4*b4, cmul(w, cmul(w2, Pb3))))

  # E = K w^(m-1), Ē = conj(K) w̄^(m-2)
  K   = (kn, ks)
  wm2 = cpow(w, m-2)
  E   = cmul(K, cmul(wm2, w))
  E_w = cscale(m-1, cmul(K, wm2))
  Eb  = cmul(cconj(K), cconj(wm2))
  if m > 2
    Eb_wb = cscale(m-2, cmul(cconj(K), cconj(cpow(w, m-3))))
  else
    Eb_wb = (zero(Eb[1]), zero(Eb[2]))
  end

  G    = cadd(cmul(E, A), cmul(Eb, B))
  G_w  = cadd(cadd(cmul(E_w, A), cmul(E, A_w)), cmul(Eb, B_w))
  G_wb = cadd(cadd(cmul(E, A_wb), cmul(Eb_wb, B)), cmul(Eb, B_wb))
  G_P  = cadd(cmul(E, A_P), cmul(Eb, B_P))
  G_Pb = cadd(cmul(E, A_Pb), cmul(Eb, B_Pb))

  # ∂W/∂x + i ∂W/∂y = P̄ G_w̄ + P conj(G_w);  ∂W/∂x' + i ∂W/∂y' = G + P̄ G_P̄ + P conj(G_P)
  gq = cadd(cmul(Pb, G_wb), cmul(P, cconj(G_w)))
  gp = cadd(cadd(G, cmul(Pb, G_Pb)), cmul(P, cconj(G_P)))
  return gq[1], gq[2], gp[1], gp[2]
end

# Longitudinal field integral Σ Im(K w^m)/m! across the edge, used for the spin rotation.
@inline function multipole_fringe_bs(mm, kn, ks, x, y)
  bs = zero(x*kn[1])
  for j in 1:length(mm)
    m = mm[j]
    if m >= 2
      fac = 1.0
      for k in 2:m
        fac *= k
      end
      Kw = cmul((kn[j], ks[j]), cpow((x, y), m))
      bs = bs + Kw[2]/fac
    end
  end
  return bs
end

@makekernel fastgtpsa=true function multipole_fringe!(i, coords::Coords, mm, kn, ks, sign)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  rel_p = 1 + v[i,PZI]

  x0  = v[i,XI]
  px0 = v[i,PXI]
  y0  = v[i,YI]
  py0 = v[i,PYI]
  x1  = x0
  px1 = px0
  y1  = y0
  py1 = py0
  dz  = zero(x0*kn[1])

  for _ in 1:MULTIPOLE_FRINGE_ITERATIONS
    xm  = (x0 + x1)/2
    pxm = (px0 + px1)/2
    ym  = (y0 + y1)/2
    pym = (py0 + py1)/2

    ps2 = rel_p*rel_p - pxm*pxm - pym*pym
    good_momenta = (ps2 > 0)
    coords.state[i] = vifelse(!good_momenta & alive, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)
    ps = sqrt(vifelse(good_momenta, ps2, one(ps2)))
    xp = pxm/ps
    yp = pym/ps

    Wx  = zero(dz)
    Wy  = zero(dz)
    Wxp = zero(dz)
    Wyp = zero(dz)
    for j in 1:length(mm)
      if mm[j] >= 2
        gx, gy, gxp, gyp = multipole_fringe_gradient(mm[j], kn[j], ks[j], xm, ym, xp, yp)
        Wx  = Wx  + gx
        Wy  = Wy  + gy
        Wxp = Wxp + gxp
        Wyp = Wyp + gyp
      end
    end

    x1  = x0  + sign*(Wxp*(1 + xp*xp) + Wyp*xp*yp)/ps
    y1  = y0  + sign*(Wxp*xp*yp + Wyp*(1 + yp*yp))/ps
    px1 = px0 - sign*Wx
    py1 = py0 - sign*Wy
    dz  = -sign*rel_p*(Wxp*xp + Wyp*yp)/(ps*ps)
  end

  v[i,XI]  = vifelse(alive, x1, v[i,XI])
  v[i,PXI] = vifelse(alive, px1, v[i,PXI])
  v[i,YI]  = vifelse(alive, y1, v[i,YI])
  v[i,PYI] = vifelse(alive, py1, v[i,PYI])
  v[i,ZI]  = vifelse(alive, v[i,ZI] + dz, v[i,ZI])
end
