@inline function thick_pure_bsolenoid(tm::Exact, kc, p_over_q_ref, bunch, bm0, L)
  Ksol, _ = get_strengths(bm0, L, p_over_q_ref)
  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)
  return push(kc, make_kernel_call(BeamTracking.exact_solenoid!, (Ksol, beta_0, gamsqr_0, tilde_m, gyromagnetic_anomaly(bunch.species), L)))
end

@inline function drift(tm::Exact, kc, p_over_q_ref, bunch, L)
  tilde_m, gamsqr_0, beta_0 = drift_params(bunch.species, p_over_q_ref)
  return push(kc, make_kernel_call(BeamTracking.exact_drift!, (0, beta_0, gamsqr_0, tilde_m, L)))
end

@inline function thick_bend_pure_bdipole(tm::Exact, kc, p_over_q_ref, bunch, bendparams, bm1, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  if tm.fringe_at == Fringe.BothEnds || tm.fringe_at == Fringe.EntranceEnd
    e1 = bendparams.e1
  else
    e1 = 0
  end
  if tm.fringe_at == Fringe.BothEnds || tm.fringe_at == Fringe.ExitEnd
    e2 = bendparams.e2
  else
    e2 = 0
  end
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  Kn0, Ks0 = get_strengths(bm1, L, p_over_q_ref)
  Ks0 == 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, _, beta_0 = drift_params(bunch.species, p_over_q_ref)
  return push(kc, make_kernel_call(BeamTracking.exact_bend_with_rotation!, (e1, e2, theta, gyromagnetic_anomaly(bunch.species), g, Kn0, w, w_inv, tilde_m, beta_0, L)))
end

@inline function thick_pure_bdipole(tm::Exact, kc, p_over_q_ref, bunch, bm1, L)
  Kn0, Ks0 = get_strengths(bm1, L, p_over_q_ref)
  Kn = sqrt(Kn0^2 + Ks0^2)
  tilt = atan2(Ks0, Kn0)
  w = rot_quaternion(0,0,tilt)
  w_inv = inv_rot_quaternion(0,0,tilt)
  tilde_m, _, beta_0 = drift_params(bunch.species, p_over_q_ref)
  return push(kc, make_kernel_call(BeamTracking.exact_bend_with_rotation!, (0, 0, 0, gyromagnetic_anomaly(bunch.species), 0, Kn, w, w_inv, tilde_m, beta_0, L)))
end

@inline function thick_bend_no_field(tm::Exact, kc, p_over_q_ref, bunch, bendparams, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = rot_quaternion(0,0,-tilt)
  w_inv = inv_rot_quaternion(0,0,-tilt)
  theta = g * L
  tilde_m, _, beta_0 = drift_params(bunch.species, p_over_q_ref)
  return push(kc, make_kernel_call(BeamTracking.exact_bend_with_rotation!, (e1, e2, theta, gyromagnetic_anomaly(bunch.species), g, 0, w, w_inv, tilde_m, beta_0, L)))
end
#=
      # EMultipole
      n_multipoles = get_n_multipoles(ele.emultipoleparams)
      if n_multipoles != 1 || emultipoleparams.order[1] != 1
        error(" EMultipoleParams is currently only defined for the cyclic order 1 electric multipole (dipole)")
      else

      end
      em1 = first(ele.emultipoleparams)
      mm = em1.order
      kn, ks = get_e_strengths(em1, L)
      kE = sqrt(kn^2 + ks^2)/(p_over_q_ref * C_LIGHT)
      tilt = atan2(ks, kn)
      if tilt == 0
        w = nothing
        w_inv = nothing
      else
        w = rot_quaternion(0, 0, tilt)
        w_inv = inv_rot_quaternion(0, 0, tilt)
      end
      q = chargeof(bunch.species)
      mc2 = massof(bunch.species)
      E_ref = mc2/tilde_m/beta_0
      a = gyromagnetic_anomaly(bunch.species)
      edge_params = (a, tilde_m, kE, w, w_inv)
      radiation_params = ele.universalparams.rad_damping ? (q, mc2, E_ref) : nothing
      params = (radiation_params, kE, beta_0, tilde_m, w, w_inv, a)
      if isprimitivetype(eltype(bunch.coords.v)) && ele.universalparams.rad_fluctuations
        photon_params = (BeamTracking.exact_elsep!, get_backend(bunch.coords.v), q, mc2, E_ref, kE, w, w_inv)
      else
        photon_params = nothing
      end
      #tm+
      kc = push(kc, make_kernel_call(BeamTracking.bmad_to_mad!, (beta_0, tilde_m, 0)))
      kc = push(kc, integration_launcher(BeamTracking.elsep_rad!, params, photon_params, tm, edge_params, L))
      kc = push(kc, make_kernel_call(BeamTracking.mad_to_bmad!, (beta_0, tilde_m, 0)))
      kc = push_transforms_out(kc, make_kernel_call(BeamTracking.callback_electric!, (beta_0, tilde_m, kE, Val{false}())))
      kc = push_transforms_in(kc, make_kernel_call(BeamTracking.callback_electric!, (beta_0, tilde_m, kE, Val{true}())))
      return kc
    end
    =#