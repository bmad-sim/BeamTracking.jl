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
  order = tm.order
  L = universalparams.L
  n_steps = universalparams.n_steps
  ds_step = L / n_steps
  fin  = fringe_in(tm.fringe_at)
  fout = fringe_out(tm.fringe_at)
  if order == 2
    return make_kernel_call(BeamTracking.order_two_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, Val{tm.opt_scheme}(), L))
  elseif order == 4
    return make_kernel_call(BeamTracking.order_four_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, Val{tm.opt_scheme}(), L))
  elseif order == 6
    return make_kernel_call(BeamTracking.order_six_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, Val{tm.opt_scheme}(), L))
  elseif order == 8
    return make_kernel_call(BeamTracking.order_eight_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, Val{tm.opt_scheme}(), L))
  elseif order == 10
    return make_kernel_call(BeamTracking.order_ten_integrator!, (ker, params, photon_params, ds_step, n_steps, edge_params, fin, fout, Val{tm.opt_scheme}(), L))
  else
    error("Integration order $order not supported")
  end
end
