# =========== HELPER FUNCTIONS ============= #
@inline function fringe_in(f::Fringe.T)
  if f == Fringe.BothEnds || f == Fringe.EntranceEnd
    return Val{true}()
  else
    return Val{false}()
  end
end

@inline function fringe_out(f::Fringe.T)
  if f == Fringe.BothEnds || f == Fringe.ExitEnd
    return Val{true}()
  else
    return Val{false}()
  end
end

@inline function integration_launcher(ker, params, photon_params, edge_params, tm, universalparams)
  if tm isa AbstractSymplectic
    order = tm.order
    opt_scheme = tm.opt_scheme ? Val{true}() : Val{false}()
  elseif tm isa Exact
    order = 2
    opt_scheme = Val{false}()
  else
    error("integration_launcher only compatible with tracking methods <: AbstractSymplectic, or Exact")
  end
  L = universalparams.L
  n_steps = universalparams.n_steps
  ds_step = universalparams.ds_step
  L_slice = universalparams.L_slice # L, or the length of a slice of the element
  @assert ds_step ≈ L_slice / n_steps "Inconsistent ds_step/n_steps: !(ds_step ≈ L_slice / n_steps) -> !($ds_step ≈ $L_slice / $n_steps) "
  fin  = fringe_in(universalparams.fringe_at)
  fout = fringe_out(universalparams.fringe_at)
  if order == 2
    return make_kernel_call(BeamTracking.order_two_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, opt_scheme, L))
  elseif order == 4
    return make_kernel_call(BeamTracking.order_four_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, opt_scheme, L))
  elseif order == 6
    return make_kernel_call(BeamTracking.order_six_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, opt_scheme, L))
  elseif order == 8
    return make_kernel_call(BeamTracking.order_eight_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, opt_scheme, L))
  elseif order == 10
    return make_kernel_call(BeamTracking.order_ten_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, opt_scheme, L))
  else
    error("Integration order $order not supported")
  end
end
