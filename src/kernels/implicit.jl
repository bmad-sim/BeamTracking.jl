@inline function implicit_spin_rad!(i, coords::Coords, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds, ::Val{radiation_damping}) where {radiation_damping}
  @inbounds begin @FastGTPSA begin
    # rotate into curvature plane

    # s += ds/2
    #if radiation_damping
      #deterministic_radiation!()
    #end

    #if isnothing(coords.q)
      #rotate_spin_implicit!()
    #end

    # go into energy
    implicit_step!(i, coords, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
    # go out of energy

    #if isnothing(coords.q)
      #rotate_spin_implicit!()
    #end

    #if radiation_damping
      #deterministic_radiation!()
    #end
    # s += ds/2

    # rotate out of curvature plane
  end end
  return nothing
end


@makekernel fastgtpsa=true function implicit_step!(i, coords::Coords, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
  v = coords.v
  v_new = (v[XI], v[PXI], v[YI], v[PYI], v[ZI], v[PZI])
  alive = (coords.state[i] == STATE_ALIVE)

  x_new = find_root_x(v_new, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
  v_new = (x_new[1], v[PXI], x_new[2], v[PYI], x_new[3], v[PZI])
  p_new = (v_new[PXI], v_new[PYI], v_new[PZI]) .- (ds .* dH_dx(v_new, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized))

  v_new = (x_new[1], p_new[1], x_new[2], p_new[2], x_new[3], p_new[3])
  p_new = find_root_p(v_new, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
  v_new = (x_new[1], p_new[1], x_new[2], p_new[2], x_new[3], p_new[3])
  x_new = (v_new[XI], v_new[YI], v_new[ZI]) .+ (ds .* dH_dp(v_new, s, beta_0, tilde_m, g, potential, p_over_q_ref, normalized))

  v[i, XI]  = vifelse(alive, v_new[XI],  v[i, XI])
  v[i, PXI] = vifelse(alive, v_new[PXI], v[i, PXI])
  v[i, YI]  = vifelse(alive, v_new[YI],  v[i, YI])
  v[i, PYI] = vifelse(alive, v_new[PYI], v[i, PYI])
  v[i, ZI]  = vifelse(alive, v_new[ZI],  v[i, ZI])
  v[i, PZI] = vifelse(alive, v_new[PZI], v[i, PZI])
  # Need an alive check at the end
end


@inline function find_root_x(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
  @inbounds begin @FastGTPSA begin
    T = eltype(v)
    ε = eps(T)
    N_max = 100
    N = 1
    conv = false
    x =  (v[XI], v[YI], v[ZI])
    x0 = (v[XI], v[YI], v[ZI])
    while !conv && N <= N_max
      v_new = (x[1], v[PXI], x[2], v[PYI], x[3], v[PZI])
      hess = mixed_hessian_H(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized)
      J = (1 - ds*hess[1],    -ds*hess[4],    -ds*hess[7],
              -ds*hess[2], 1 - ds*hess[5],    -ds*hess[8],
              -ds*hess[3],    -ds*hess[6], 1 - ds*hess[9])
      F = x .- x0 .- (ds .* dH_dp(v_new, s, beta_0, tilde_m, g, potential, p_over_q_ref, normalized))
      sol = solve_3x3_cramer(J, -1 .* F)
      conv = all((sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3])/(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]) < ε)
      x = x .+ sol
      N += 1
    end
    if !conv && N == N_max
      @warn "Implicit integrator's Newton search did not converge in $N_max iterations"
      end
    return x
  end end
end


@inline function find_root_p(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized, ds)
  @inbounds begin @FastGTPSA begin
    T = eltype(v)
    ε = eps(T)
    N_max = 100
    N = 1
    conv = false
    p =  (v[PXI], v[PYI], v[PZI])
    p0 = (v[PXI], v[PYI], v[PZI])
    while !conv && N <= N_max
      v_new = (v[XI], p[1], v[YI], p[2], v[ZI], p[3])
      hess = mixed_hessian_H(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized)
      J = (1 - ds*hess[1],    -ds*hess[2],    -ds*hess[3],
              -ds*hess[4], 1 - ds*hess[5],    -ds*hess[6],
              -ds*hess[7],    -ds*hess[8], 1 - ds*hess[9])
      F = p .- p0 .+ (ds .* dH_dx(v_new, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, normalized))
      sol = solve_3x3_cramer(J, -1 .* F)
      conv =  all((sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3])/(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]) < ε)
      p = p .+ sol
      N =+ 1
    end
    if !conv && N == N_max
      @warn "Implicit integrator's Newton search did not converge in $N_max iterations"
      end
    return p
  end end
end


"""
Returns the position derivatives of the Hamiltonian.
"""
@inline function dH_dx(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, ::Val{normalized}) where {normalized}
  @inbounds begin @FastGTPSA begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    phi, ax, ay, az = potential(v[XI], v[YI], s, t) # as is called az because as is a keyword
    derivatives = jac(v[XI], v[YI], s, t)
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
      az =  az/p_over_q_ref

      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/p_over_q_ref/C_LIGHT, derivatives[5]/p_over_q_ref, derivatives[9]/p_over_q_ref, derivatives[13]/p_over_q_ref
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/p_over_q_ref/C_LIGHT, derivatives[6]/p_over_q_ref, derivatives[10]/p_over_q_ref, derivatives[14]/p_over_q_ref
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3]/p_over_q_ref/C_LIGHT, derivatives[7]/p_over_q_ref, derivatives[11]/p_over_q_ref, derivatives[15]/p_over_q_ref
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/p_over_q_ref/C_LIGHT, derivatives[8]/p_over_q_ref, derivatives[12]/p_over_q_ref, derivatives[16]/p_over_q_ref
    else
      phi = phi/C_LIGHT
      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/C_LIGHT, derivatives[5], derivatives[9], derivatives[13]
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/C_LIGHT, derivatives[6], derivatives[10], derivatives[14]
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3]/C_LIGHT, derivatives[7], derivatives[11], derivatives[15]
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/C_LIGHT, derivatives[8], derivatives[12], derivatives[16]
    end
    dphi_dz = beta_0*dphi_ds - dphi_dt/C_LIGHT
    dax_dz =  beta_0*dax_ds  - dax_dt/C_LIGHT
    day_dz =  beta_0*day_ds  - day_dt/C_LIGHT
    daz_dz =  beta_0*daz_ds  - daz_dt/C_LIGHT

    px = v[PXI] - ax
    py = v[PYI] - ay
    rel_p = v[PZI] + 1/beta_0 - phi
    
    ps2 = rel_p*rel_p - tilde_m*tilde_m - px*px - py*py
    good_momenta = (ps2 > 0)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))

    dH_dx = h*(rel_p*dphi_dx - px*dax_dx - py*day_dx)/ps - daz_dx - g*ps
    dH_dy = h*(rel_p*dphi_dy - px*dax_dy - py*day_dy)/ps - daz_dy
    dH_dz = h*(rel_p*dphi_dz - px*dax_dz - py*day_dz)/ps - daz_dz

    return (dH_dx, dH_dy, dH_dz)
  end end
end


"""
Returns the momentum derivatives of the Hamiltonian.
"""
@inline function dH_dp(v, s, beta_0, tilde_m, g, potential, p_over_q_ref, ::Val{normalized}) where {normalized}
  @inbounds begin @FastGTPSA begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    phi, ax, ay, az = potential(v[XI], v[YI], s, t) # as is called az because as is a keyword
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
      az =  az/p_over_q_ref
    else
      phi = phi/C_LIGHT
    end
    px = v[PXI] - ax
    py = v[PYI] - ay
    rel_p = v[PZI] + 1/beta_0 - phi
    
    ps2 = rel_p*rel_p - tilde_m*tilde_m - px*px - py*py
    good_momenta = (ps2 > 0)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))
    dH_dpx =  h*px/ps
    dH_dpy =  h*py/ps
    dH_dpz = -h*rel_p/ps + 1/beta_0
    
    return (dH_dpx, dH_dpy, dH_dpz)
  end end
end


"""
Returns the mixed position-momentum second derivatives of the Hamiltonian.
"""
@inline function mixed_hessian_H(v, s, beta_0, tilde_m, g, potential, jac, p_over_q_ref, ::Val{normalized}) where {normalized}
  @inbounds begin @FastGTPSA begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    phi, ax, ay, az = potential(v[XI], v[YI], s, t) # as is called az because as is a keyword
    derivatives = jac(v[XI], v[YI], s, t)
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
      az =  az/p_over_q_ref

      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/p_over_q_ref/C_LIGHT, derivatives[5]/p_over_q_ref, derivatives[9]/p_over_q_ref, derivatives[13]/p_over_q_ref
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/p_over_q_ref/C_LIGHT, derivatives[6]/p_over_q_ref, derivatives[10]/p_over_q_ref, derivatives[14]/p_over_q_ref
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3]/p_over_q_ref/C_LIGHT, derivatives[7]/p_over_q_ref, derivatives[11]/p_over_q_ref, derivatives[15]/p_over_q_ref
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/p_over_q_ref/C_LIGHT, derivatives[8]/p_over_q_ref, derivatives[12]/p_over_q_ref, derivatives[16]/p_over_q_ref
    else
      phi = phi/C_LIGHT
      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/C_LIGHT, derivatives[5], derivatives[9], derivatives[13]
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/C_LIGHT, derivatives[6], derivatives[10], derivatives[14]
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3]/C_LIGHT, derivatives[7], derivatives[11], derivatives[15]
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/C_LIGHT, derivatives[8], derivatives[12], derivatives[16]
    end
    dphi_dz = beta_0*dphi_ds - dphi_dt/C_LIGHT
    dax_dz =  beta_0*dax_ds -  dax_dt/C_LIGHT
    day_dz =  beta_0*day_ds -  day_dt/C_LIGHT
    daz_dz =  beta_0*daz_ds -  daz_dt/C_LIGHT

    px = v[PXI] - ax
    py = v[PYI] - ay
    rel_p = v[PZI] + 1/beta_0 - phi
    
    ps2 = rel_p*rel_p - tilde_m*tilde_m - px*px - py*py
    good_momenta = (ps2 > 0)
    ps2_1 = one(ps2)
    ps = sqrt(vifelse(good_momenta, ps2, ps2_1))
    h_over_ps3 = h/(ps2*ps)

    middle_factor_x = h_over_ps3*(rel_p*dphi_dx - px*dax_dx - py*day_dx)
    middle_factor_y = h_over_ps3*(rel_p*dphi_dy - px*dax_dy - py*day_dy)
    middle_factor_z = h_over_ps3*(rel_p*dphi_dz - px*dax_dz - py*day_dz)

    d2H_dxdpx = px*middle_factor_x - h*dax_dx/ps + g*px/ps
    d2H_dxdpy = py*middle_factor_x - h*day_dx/ps + g*py/ps
    d2H_dxdpz = -rel_p*middle_factor_x + h*dphi_dx/ps - g*rel_p/ps

    d2H_dydpx = px*middle_factor_y - h*dax_dy/ps
    d2H_dydpy = py*middle_factor_y - h*day_dy/ps
    d2H_dydpz = -rel_p*middle_factor_y + h*dphi_dy/ps

    d2H_dzdpx = px*middle_factor_z - h*dax_dz/ps
    d2H_dzdpy = py*middle_factor_z - h*day_dz/ps
    d2H_dzdpz = -rel_p*middle_factor_z + h*dphi_dz/ps

    return (d2H_dxdpx, d2H_dxdpy, d2H_dxdpz, d2H_dydpx, d2H_dydpy, d2H_dydpz, d2H_dzdpx, d2H_dzdpy, d2H_dzdpz)
  end end
end


"""
Solves Ax = y for x using Cramer's rule, where A is a 3x3 matrix and y is a 3-vector.
"""
@inline function solve_3x3_cramer(A, y)
  @inbounds begin
    y1, y2, y3 = y
    a, b, c, d, e, f, g, h, i = A

    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    d1 = y1*(e*i - f*h) - b*(y2*i - f*y3) + c*(y2*h - e*y3)
    d2 = a*(y2*i - f*y3) - y1*(d*i - f*g) + c*(d*y3 - y2*g)
    d3 = a*(e*y3 - y2*h) - b*(d*y3 - y2*g) + y1*(d*h - e*g)

    return (d1/det, d2/det, d3/det)
  end
end