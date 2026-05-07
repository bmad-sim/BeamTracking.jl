function implicit_integrator!(i, coords::Coords, s, radiation_params, beta_0, tilde_m, a, g, w, w_inv, potential_and_jac::U, potential_params, p_over_q_ref, normalized, L) where {U}
  @inbounds begin
    rotation!(i, coords, w, 0)

    s += L / 2

    if !isnothing(radiation_params)
      q, mc2, E_ref = radiation_params
      deterministic_radiation_implicit!(i, coords, s, q, mc2, E_ref, g, potential_and_jac, potential_params, p_over_q_ref, normalized, L / 2)
    end

    if !isnothing(coords.q)
      rotate_spin_implicit!(i, coords, s, a, g, beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, L / 2)
    end

    implicit_step!(i, coords, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, normalized, L)

    if !isnothing(coords.q)
      rotate_spin_implicit!(i, coords, s, a, g, beta_0, tilde_m, potential_and_jac, potential_params, p_over_q_ref, normalized, L / 2)
    end

    if !isnothing(radiation_params)
      deterministic_radiation_implicit!(i, coords, s, q, mc2, E_ref, g, potential_and_jac, potential_params, p_over_q_ref, normalized, L / 2)
    end

    rotation!(i, coords, w_inv, 0)
  end
  return nothing
end


function implicit_step!(i, coords::Coords, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, ::Val{normalized}, ds) where {U, normalized}
  @inbounds begin
    v = coords.v
    alive_at_start = (coords.state[i] == STATE_ALIVE)

    v_orig = scalar.((v[i,XI], v[i,PXI], v[i,YI], v[i,PYI], v[i,ZI], v[i,PZI]))
    v_new = v_orig

    x_new = scalar.(find_root_x(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}(), ds/2))
    v_new = (x_new[1], v_new[PXI], x_new[2], v_new[PYI], x_new[3], v_new[PZI])

    p_new = (v_new[PXI], v_new[PYI], v_new[PZI]) .- (ds/2 .* scalar.(dH_dx(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())))
    v_new = (v_new[XI], p_new[1], v_new[YI], p_new[2], v_new[ZI], p_new[3])

    if eltype(v) <: TPS
      nn = TPSAInterface.ndiffs(v[i,1])
      f = similar(eltype(v), nn)
      for j in 1:nn
        f[j] = zero(v[i,XI])
        TPSAInterface.seti!(f[j], 1, j)
      end
  
      TPSAInterface.seti!(f[XI],   v_new[XI],   0)
      TPSAInterface.seti!(f[YI],   v_new[YI],   0)
      TPSAInterface.seti!(f[ZI],   v_new[ZI],   0)
      TPSAInterface.seti!(f[PXI],  v_orig[PXI], 0)
      TPSAInterface.seti!(f[PYI],  v_orig[PYI], 0)
      TPSAInterface.seti!(f[PZI],  v_orig[PZI], 0)

      d1 = ds/2 .* dH_dx(f, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())
      d2 = ds/2 .* dH_dp(f, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())
    
      f[XI]  = f[XI]  - v_orig[XI] - d2[1]
      f[YI]  = f[YI]  - v_orig[YI] - d2[2]
      f[ZI]  = f[ZI]  - v_orig[ZI] - d2[3]
      f[PXI] = f[PXI] - v_new[PXI] - d1[1]
      f[PYI] = f[PYI] - v_new[PYI] - d1[2]
      f[PZI] = f[PZI] - v_new[PZI] - d1[3]

      f2 = zero(f)
      inds = zeros(Cint, nn)
      inds[XI] = 1
      inds[YI] = 1
      inds[ZI] = 1
      # Long term solution would be to move pminv into TPSAInterface,
      # For now use GTPSA directly
      GTPSA.pminv!(nn, f, 6, f2, inds)

      f3 = zero(f)
      for j in 1:6
        f3[j] = TPSAInterface.cutord(v[i,j], 0)
      end
      for j in 7:nn
        f3[j] = 0
        TPSAInterface.seti!(f3[j], 1, j)
      end
      # Uses GTPSA's compose ∘, long term solution would be to move 
      # to TPSAInterface
      v_final = Tuple(v_new .+ (f2 ∘ f3)[1:6])
      v_orig = v_new
    elseif eltype(v) <: ForwardDiff.Dual
      A, B, C = let vx = v_new[XI],
                    vy = v_new[YI],
                    vz = v_new[ZI],
                    vpx = v_orig[PXI],
                    vpy = v_orig[PYI],
                    vpz = v_orig[PZI]
        hess = ds/2 .* mixed_hessian_H((vx, vpx, vy, vpy, vz, vpz), s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())
        A = SA[hess[1] hess[2] hess[3];
               hess[4] hess[5] hess[6];
               hess[7] hess[8] hess[9]]
        B = ds/2 .* ForwardDiff.jacobian(p -> SVector{3}(dH_dp(SA[vx, p[1], vy, p[2], vz, p[3]], s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())), 
        SA[vpx, vpy, vpz])
        C = ds/2 .* ForwardDiff.jacobian(x -> SVector{3}(dH_dx(SA[x[1], vpx, x[2], vpy, x[3], vpz], s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())), 
        SA[vx, vy, vz])
        A, B, C
      end
      id = SA[1 0 0;
              0 1 0;
              0 0 1]
      Mxx = inv(id .- A')
      Mxp = Mxx * B
      Mpx = -C * Mxx
      Mpp = I - A - (C * Mxx * B)
      jac = SA[Mxx[1,1] Mxp[1,1] Mxx[1,2] Mxp[1,2] Mxx[1,3] Mxp[1,3];
              Mpx[1,1] Mpp[1,1] Mpx[1,2] Mpp[1,2] Mpx[1,3] Mpp[1,3];
              Mxx[2,1] Mxp[2,1] Mxx[2,2] Mxp[2,2] Mxx[2,3] Mxp[2,3];
              Mpx[2,1] Mpp[2,1] Mpx[2,2] Mpp[2,2] Mpx[2,3] Mpp[2,3];
              Mxx[3,1] Mxp[3,1] Mxx[3,2] Mxp[3,2] Mxx[3,3] Mxp[3,3];
              Mpx[3,1] Mpp[3,1] Mpx[3,2] Mpp[3,2] Mpx[3,3] Mpp[3,3]]
      r1 = ForwardDiff.partials(v[i,1]).values
      r2 = ForwardDiff.partials(v[i,2]).values
      r3 = ForwardDiff.partials(v[i,3]).values
      r4 = ForwardDiff.partials(v[i,4]).values
      r5 = ForwardDiff.partials(v[i,5]).values
      r6 = ForwardDiff.partials(v[i,6]).values
      jac_orig = SA[r1[1] r1[2] r1[3] r1[4] r1[5] r1[6];
                    r2[1] r2[2] r2[3] r2[4] r2[5] r2[6];
                    r3[1] r3[2] r3[3] r3[4] r3[5] r3[6];
                    r4[1] r4[2] r4[3] r4[4] r4[5] r4[6];
                    r5[1] r5[2] r5[3] r5[4] r5[5] r5[6];
                    r6[1] r6[2] r6[3] r6[4] r6[5] r6[6]]
      jac_new = jac * jac_orig
      T = eltype(v).parameters[1]
      new_x  = ForwardDiff.Dual{T}(v_new[XI],  Tuple(jac_new[XI,:]))
      new_y  = ForwardDiff.Dual{T}(v_new[YI],  Tuple(jac_new[YI,:]))
      new_z  = ForwardDiff.Dual{T}(v_new[ZI],  Tuple(jac_new[ZI,:]))
      new_px = ForwardDiff.Dual{T}(v_new[PXI], Tuple(jac_new[PXI,:]))
      new_py = ForwardDiff.Dual{T}(v_new[PYI], Tuple(jac_new[PYI,:]))
      new_pz = ForwardDiff.Dual{T}(v_new[PZI], Tuple(jac_new[PZI,:]))
      v_final = (new_x, new_px, new_y, new_py, new_z, new_pz)
      v_orig = v_new
    else
      v_final = v_new
    end

    p_new = scalar.(find_root_p(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}(), ds/2))
    v_new = (v_new[XI], p_new[1], v_new[YI], p_new[2], v_new[ZI], p_new[3])

    x_new = (v_new[XI], v_new[YI], v_new[ZI]) .+ (ds/2 .* scalar.(dH_dp(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())))
    v_new = (x_new[1], v_new[PXI], x_new[2], v_new[PYI], x_new[3], v_new[PZI])

    if eltype(v) <: TPS
      for j in 1:nn
        TPSAInterface.clear!(f[j])
        TPSAInterface.seti!(f[j], 1, j)
      end
      TPSAInterface.seti!(f[XI],  v_orig[XI], 0)
      TPSAInterface.seti!(f[YI],  v_orig[YI], 0)
      TPSAInterface.seti!(f[ZI],  v_orig[ZI], 0)
      TPSAInterface.seti!(f[PXI], v_new[PXI], 0)
      TPSAInterface.seti!(f[PYI], v_new[PYI], 0)
      TPSAInterface.seti!(f[PZI], v_new[PZI], 0)
      d1 = ds/2 .* dH_dx(f, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())
      d2 = ds/2 .* dH_dp(f, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())

      f[XI]  = f[XI]  - v_new[XI]   + d2[1]
      f[YI]  = f[YI]  - v_new[YI]   + d2[2]
      f[ZI]  = f[ZI]  - v_new[ZI]   + d2[3]
      f[PXI] = f[PXI] - v_orig[PXI] + d1[1]
      f[PYI] = f[PYI] - v_orig[PYI] + d1[2]
      f[PZI] = f[PZI] - v_orig[PZI] + d1[3]

      inds[XI]  = 0
      inds[PXI] = 1
      inds[YI]  = 0
      inds[PYI] = 1
      inds[ZI]  = 0
      inds[PZI] = 1
      GTPSA.pminv!(nn, f, 6, f2, inds)

      for j in 1:6
        f3[j] = TPSAInterface.cutord(v_final[j], 0)
      end
      for j in 7:nn
        TPSAInterface.clear!(f3[j])
        TPSAInterface.seti!(f3[j], 1, j)
      end
      v_final = Tuple(v_new .+ (f2 ∘ f3)[1:6])
    elseif eltype(v) <: ForwardDiff.Dual
      A, B, C = let vx = v_orig[XI],
                    vy = v_orig[YI],
                    vz = v_orig[ZI],
                    vpx = v_new[PXI],
                    vpy = v_new[PYI],
                    vpz = v_new[PZI]
        hess = ds/2 .* mixed_hessian_H((vx, vpx, vy, vpy, vz, vpz), s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())
        A = SA[hess[1] hess[2] hess[3];
               hess[4] hess[5] hess[6];
               hess[7] hess[8] hess[9]]
        B = ds/2 .* ForwardDiff.jacobian(p -> SVector{3}(dH_dp(SA[vx, p[1], vy, p[2], vz, p[3]], s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())), 
        SA[vpx, vpy, vpz])
        C = ds/2 .* ForwardDiff.jacobian(x -> SVector{3}(dH_dx(SA[x[1], vpx, x[2], vpy, x[3], vpz], s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, Val{normalized}())), 
        SA[vx, vy, vz])
        A, B, C
      end
      Mpp = inv(id .+ A)
      Mxp = B * Mpp
      Mpx = -Mpp * C
      Mxx = id .+ A' .- (B * Mpp * C)
      jac = SA[Mxx[1,1] Mxp[1,1] Mxx[1,2] Mxp[1,2] Mxx[1,3] Mxp[1,3];
               Mpx[1,1] Mpp[1,1] Mpx[1,2] Mpp[1,2] Mpx[1,3] Mpp[1,3];
               Mxx[2,1] Mxp[2,1] Mxx[2,2] Mxp[2,2] Mxx[2,3] Mxp[2,3];
               Mpx[2,1] Mpp[2,1] Mpx[2,2] Mpp[2,2] Mpx[2,3] Mpp[2,3];
               Mxx[3,1] Mxp[3,1] Mxx[3,2] Mxp[3,2] Mxx[3,3] Mxp[3,3];
               Mpx[3,1] Mpp[3,1] Mpx[3,2] Mpp[3,2] Mpx[3,3] Mpp[3,3]]
      r1 = ForwardDiff.partials(v_final[1]).values
      r2 = ForwardDiff.partials(v_final[2]).values
      r3 = ForwardDiff.partials(v_final[3]).values
      r4 = ForwardDiff.partials(v_final[4]).values
      r5 = ForwardDiff.partials(v_final[5]).values
      r6 = ForwardDiff.partials(v_final[6]).values
      jac_orig = SA[r1[1] r1[2] r1[3] r1[4] r1[5] r1[6];
                    r2[1] r2[2] r2[3] r2[4] r2[5] r2[6];
                    r3[1] r3[2] r3[3] r3[4] r3[5] r3[6];
                    r4[1] r4[2] r4[3] r4[4] r4[5] r4[6];
                    r5[1] r5[2] r5[3] r5[4] r5[5] r5[6];
                    r6[1] r6[2] r6[3] r6[4] r6[5] r6[6]]
      jac_new = jac * jac_orig
      new_x  = ForwardDiff.Dual{T}(v_new[XI],  Tuple(jac_new[XI,:]))
      new_y  = ForwardDiff.Dual{T}(v_new[YI],  Tuple(jac_new[YI,:]))
      new_z  = ForwardDiff.Dual{T}(v_new[ZI],  Tuple(jac_new[ZI,:]))
      new_px = ForwardDiff.Dual{T}(v_new[PXI], Tuple(jac_new[PXI,:]))
      new_py = ForwardDiff.Dual{T}(v_new[PYI], Tuple(jac_new[PYI,:]))
      new_pz = ForwardDiff.Dual{T}(v_new[PZI], Tuple(jac_new[PZI,:]))
      v_final = (new_x, new_px, new_y, new_py, new_z, new_pz)
    else
      v_final = v_new
    end

    t = (s/beta_0 - v_final[ZI])/C_LIGHT
    potential, _ = potential_and_jac(v_final[XI], v_final[YI], s, t, potential_params)
    phi, ax, ay, _ = potential
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
    else
      phi = phi/C_LIGHT
    end
    rel_p = v_final[PZI] + 1/beta_0 - phi
    px = v_final[PXI] - ax
    py = v_final[PYI] - ay

    ps2 = rel_p*rel_p - tilde_m*tilde_m - px*px - py*py
    good_momenta = (ps2 > 0)
    coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
    alive = (coords.state[i] == STATE_ALIVE)

    v[i,XI]  = vifelse(alive, v_final[XI],  v[i,XI])
    v[i,PXI] = vifelse(alive, v_final[PXI], v[i,PXI])
    v[i,YI]  = vifelse(alive, v_final[YI],  v[i,YI])
    v[i,PYI] = vifelse(alive, v_final[PYI], v[i,PYI])
    v[i,ZI]  = vifelse(alive, v_final[ZI],  v[i,ZI])
    v[i,PZI] = vifelse(alive, v_final[PZI], v[i,PZI])
  end
  return nothing
end


function find_root_x(v, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, normalized, ds) where {U}
  @inbounds begin
    ε = my_eps(v[1])
    N_max = 100
    N = 1
    conv = false
    x  = (v[XI], v[YI], v[ZI])
    x0 = (v[XI], v[YI], v[ZI])
    norm_x = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
    while !conv && N <= N_max
      v_new = (x[1], v[PXI], x[2], v[PYI], x[3], v[PZI])
      hess = scalar.(mixed_hessian_H(v, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, normalized))
      J = (1 - ds*hess[1],    -ds*hess[4],    -ds*hess[7],
              -ds*hess[2], 1 - ds*hess[5],    -ds*hess[8],
              -ds*hess[3],    -ds*hess[6], 1 - ds*hess[9])
      F = x .- x0 .- (ds .* scalar.(dH_dp(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, normalized)))
      norm_F = sqrt(F[1]*F[1] + F[2]*F[2] + F[3]*F[3])
      sol = solve_3x3_cramer(J, -1 .* F)
      norm_sol = sqrt(sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3])
      conv = (all(norm_sol < ε*norm_x) || all(norm_F < ε))
      x = x .+ sol
      N += 1
    end
    if !conv
      @warn "Implicit integrator's Newton search did not converge in $N_max iterations"
    end
    return x
  end
end


function find_root_p(v, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, normalized, ds) where {U}
  @inbounds begin
    ε = my_eps(v[1])
    N_max = 100
    N = 1
    conv = false
    p  = (v[PXI], v[PYI], v[PZI])
    p0 = (v[PXI], v[PYI], v[PZI])
    norm_p = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3])
    while !conv && N <= N_max
      v_new = (v[XI], p[1], v[YI], p[2], v[ZI], p[3])
      hess = scalar.(mixed_hessian_H(v, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, normalized))
      J = (1 + ds*hess[1],     ds*hess[2],     ds*hess[3],
               ds*hess[4], 1 + ds*hess[5],     ds*hess[6],
               ds*hess[7],     ds*hess[8], 1 + ds*hess[9])
      F = p .- p0 .+ (ds .* scalar.(dH_dx(v_new, s, beta_0, tilde_m, g, potential_and_jac, potential_params, p_over_q_ref, normalized)))
      norm_F = sqrt(F[1]*F[1] + F[2]*F[2] + F[3]*F[3])
      sol = solve_3x3_cramer(J, -1 .* F)
      norm_sol = sqrt(sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3])
      conv = (all(norm_sol < ε*norm_p) || all(norm_F < ε))
      p = p .+ sol
      N += 1
    end
    if !conv
      @warn "Implicit integrator's Newton search did not converge in $N_max iterations"
    end
    return p
  end
end


"""
Returns the position derivatives of the Hamiltonian.
"""
function dH_dx(v, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, ::Val{normalized}) where {U, normalized}
  @inbounds begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    potential, derivatives = potential_and_jac(v[XI], v[YI], s, t, potential_params)
    phi, ax, ay, az = potential
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
      az =  az/p_over_q_ref

      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/p_over_q_ref/C_LIGHT, derivatives[5]/p_over_q_ref, derivatives[9]/p_over_q_ref,  derivatives[13]/p_over_q_ref
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/p_over_q_ref/C_LIGHT, derivatives[6]/p_over_q_ref, derivatives[10]/p_over_q_ref, derivatives[14]/p_over_q_ref
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/p_over_q_ref/C_LIGHT, derivatives[8]/p_over_q_ref, derivatives[12]/p_over_q_ref, derivatives[16]/p_over_q_ref
    else
      phi = phi/C_LIGHT
      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/C_LIGHT, derivatives[5], derivatives[9],  derivatives[13]
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/C_LIGHT, derivatives[6], derivatives[10], derivatives[14]
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/C_LIGHT, derivatives[8], derivatives[12], derivatives[16]
    end
    dphi_dz = -dphi_dt/C_LIGHT
    dax_dz  = -dax_dt/C_LIGHT
    day_dz  = -day_dt/C_LIGHT
    daz_dz  = -daz_dt/C_LIGHT

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
  end
end


"""
Returns the momentum derivatives of the Hamiltonian.
"""
function dH_dp(v, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, ::Val{normalized}) where {U, normalized}
  @inbounds begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    potential, _ = potential_and_jac(v[XI], v[YI], s, t, potential_params)
    phi, ax, ay, az = potential
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
  end
end


"""
Returns the mixed position-momentum second derivatives of the Hamiltonian.
"""
function mixed_hessian_H(v, s, beta_0, tilde_m, g, potential_and_jac::U, potential_params, p_over_q_ref, ::Val{normalized}) where {U, normalized}
  @inbounds begin
    h = 1 + g*v[XI]
    t = (s/beta_0 - v[ZI])/C_LIGHT

    potential, derivatives = potential_and_jac(v[XI], v[YI], s, t, potential_params)
    phi, ax, ay, az = potential
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax =  ax/p_over_q_ref
      ay =  ay/p_over_q_ref
      az =  az/p_over_q_ref

      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/p_over_q_ref/C_LIGHT, derivatives[5]/p_over_q_ref, derivatives[9]/p_over_q_ref, derivatives[13]/p_over_q_ref
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/p_over_q_ref/C_LIGHT, derivatives[6]/p_over_q_ref, derivatives[10]/p_over_q_ref, derivatives[14]/p_over_q_ref
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/p_over_q_ref/C_LIGHT, derivatives[8]/p_over_q_ref, derivatives[12]/p_over_q_ref, derivatives[16]/p_over_q_ref
    else
      phi = phi/C_LIGHT
      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/C_LIGHT, derivatives[5], derivatives[9], derivatives[13]
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/C_LIGHT, derivatives[6], derivatives[10], derivatives[14]
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/C_LIGHT, derivatives[8], derivatives[12], derivatives[16]
    end
    dphi_dz = -dphi_dt/C_LIGHT
    dax_dz =  -dax_dt/C_LIGHT
    day_dz =  -day_dt/C_LIGHT
    daz_dz =  -daz_dt/C_LIGHT

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

    return (d2H_dxdpx, d2H_dxdpy, d2H_dxdpz, 
            d2H_dydpx, d2H_dydpy, d2H_dydpz, 
            d2H_dzdpx, d2H_dzdpy, d2H_dzdpz)
  end
end


"""
Solves Ax = y for x using Cramer's rule, where A is a 3x3 matrix and y is a 3-vector.
"""
function solve_3x3_cramer(A, y)
  @inbounds begin
    y1, y2, y3 = y
    a, b, c, d, e, f, g, h, i = A

    det =  a*(e*i  - f*h)  -  b*(d*i  - f*g)  +  c*(d*h  - e*g)
    d1  = y1*(e*i  - f*h)  -  b*(y2*i - f*y3) +  c*(y2*h - e*y3)
    d2  =  a*(y2*i - f*y3) - y1*(d*i  - f*g)  +  c*(d*y3 - y2*g)
    d3  =  a*(e*y3 - y2*h) -  b*(d*y3 - y2*g) + y1*(d*h  - e*g)

    return (d1/det, d2/det, d3/det)
  end
end


"""
Returns normalized transverse vector potential and electromagnetic fields for 
implicit integrators.
"""
function implicit_fields(x, y, s, t, g, potential_and_jac::U, potential_params, p_over_q_ref, ::Val{normalized}) where {U, normalized}
  @inbounds begin @FastGTPSA begin
    h = 1 + g*x
    potential, derivatives = potential_and_jac(x, y, s, t, potential_params)
    phi, ax, ay, _ = potential
    if !normalized
      phi = phi/p_over_q_ref/C_LIGHT
      ax = ax/p_over_q_ref
      ay = ay/p_over_q_ref

      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1]/p_over_q_ref, derivatives[5]/p_over_q_ref, derivatives[9]/p_over_q_ref,  derivatives[13]/p_over_q_ref
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2]/p_over_q_ref, derivatives[6]/p_over_q_ref, derivatives[10]/p_over_q_ref, derivatives[14]/p_over_q_ref
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3]/p_over_q_ref, derivatives[7]/p_over_q_ref, derivatives[11]/p_over_q_ref, derivatives[15]/p_over_q_ref
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4]/p_over_q_ref, derivatives[8]/p_over_q_ref, derivatives[12]/p_over_q_ref, derivatives[16]/p_over_q_ref
    else
      phi = phi/C_LIGHT 
      dphi_dx, dax_dx, day_dx, daz_dx = derivatives[1], derivatives[5], derivatives[9],  derivatives[13]
      dphi_dy, dax_dy, day_dy, daz_dy = derivatives[2], derivatives[6], derivatives[10], derivatives[14]
      dphi_ds, dax_ds, day_ds, daz_ds = derivatives[3], derivatives[7], derivatives[11], derivatives[15]
      dphi_dt, dax_dt, day_dt, daz_dt = derivatives[4], derivatives[8], derivatives[12], derivatives[16]
    end
  
    # Both E and B are normalized by p_over_q_ref
    ex =  -dphi_dx - dax_dt
    ey =  -dphi_dy - day_dt
    ez = (-dphi_ds - daz_dt)/h

    bx = (daz_dy - day_ds)/h
    by = (dax_ds - daz_dx)/h
    bz =  day_dx - dax_dy

    return phi, ax, ay, ex, ey, ez, bx, by, bz
  end end
end


"""
Rotates spin for implicit integrators.
"""
function rotate_spin_implicit!(i, coords::Coords, s, a, g, beta_0, tilde_m, potential_and_jac::U, potential_params, p_over_q_ref, normalized, L) where {U}
  @inbounds begin @FastGTPSA begin
    v = coords.v

    t = (s/beta_0 - v[i,ZI])/C_LIGHT

    phi, ax, ay, ex, ey, ez, bx, by, bz = implicit_fields(v[i,XI], v[i,YI], s, t, g, potential_and_jac, potential_params, p_over_q_ref, normalized)
    e_vec = (ex, ey, ez)
    b_vec = (bx, by, bz)

    mad_to_bmad!(i, coords, beta_0, tilde_m, phi)
    rotate_spin_field!(i, coords, a, g, tilde_m, ax, ay, e_vec, b_vec, L)
    bmad_to_mad!(i, coords, beta_0, tilde_m, phi)
  end end
  return nothing
end


""" 
Applies radiation damping kick for implicit integrators.
"""
function deterministic_radiation_implicit!(i, coords::Coords, s, q, mc2, E_ref, g, potential_and_jac::U, potential_params, p_over_q_ref, normalized, L) where {U}
  @inbounds begin @FastGTPSA begin
    v = coords.v
    t = (s - v[i,ZI])/C_LIGHT # radiation is only accurate when beta_0 is approximately 1
    tilde_m = mc2/E_ref

    phi, ax, ay, ex, ey, ez, bx, by, bz = implicit_fields(v[i,XI], v[i,YI], s, t, g, potential_and_jac, potential_params, p_over_q_ref, normalized)
    e_vec = (ex, ey, ez)
    b_vec = (bx, by, bz)

    mad_to_bmad!(i, coords, 1, tilde_m, phi)
    deterministic_radiation_field!(i, coords, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)
    bmad_to_mad!(i, coords, 1, tilde_m, phi)
  end end
  return nothing
end


"""
Applies radiation diffusion kick for implicit integrators.
"""
function stochastic_radiation!(i, coords::Coords, s, ::typeof(implicit_integrator!), backend, q, mc2, E_ref, g, potential_and_jac::U, potential_params, p_over_q_ref, normalized, L) where {U}
  @inbounds begin
    v = coords.v
    t = (s - v[i,ZI])/C_LIGHT # radiation is only accurate when beta_0 is approximately 1
    tilde_m = mc2/E_ref

    phi, ax, ay, ex, ey, ez, bx, by, bz = implicit_fields(v[i,XI], v[i,YI], s, t, g, potential_and_jac, potential_params, p_over_q_ref, normalized)
    e_vec = (ex, ey, ez)
    b_vec = (bx, by, bz)

    mad_to_bmad!(i, coords, 1, tilde_m, phi)
    stochastic_radiation_field!(i, coords, backend, q, mc2, E_ref, g, ax, ay, e_vec, b_vec, L)
    bmad_to_mad!(i, coords, 1, tilde_m, phi)
  end
  return nothing
end


scalar(x::TPS) = TPSAInterface.scalar(x)
scalar(x::ForwardDiff.Dual) = ForwardDiff.value(x)
scalar(x) = x

my_eps(::SIMD.Vec{N,T}) where {N,T} = eps(T)
my_eps(::T) where {T} = eps(T)