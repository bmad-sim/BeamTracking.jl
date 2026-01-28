#
# ===============  I N T E G R A T O R S  ===============
#

@makekernel fastgtpsa=true function order_two_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, num_steps, edge_params, L)
  if !isnothing(edge_params)
    a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  for step in 1:num_steps
    ker(i, coords, params..., ds_step)
    if !isnothing(photon_params) && (step < num_steps)
      stochastic_radiation!(i, coords, photon_params..., ds_step)
    end
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  if !isnothing(edge_params)
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
  end
end


@makekernel fastgtpsa=true function order_four_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, num_steps, edge_params, L)
  w0 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  w1 =  1.3512071919596577718181151794851757586002349853515625*ds_step
  if !isnothing(edge_params)
    a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  for step in 1:num_steps
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1)
    if !isnothing(photon_params) && (step < num_steps)
      stochastic_radiation!(i, coords, photon_params..., ds_step)
    end
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  if !isnothing(edge_params)
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
  end
end


@makekernel fastgtpsa=true function order_six_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, num_steps, edge_params, L)
  w0 =  1.315186320683911169737712043570355*ds_step
  w1 = -1.17767998417887100694641568096432*ds_step
  w2 =  0.235573213359358133684793182978535*ds_step
  w3 =  0.784513610477557263819497633866351*ds_step
  if !isnothing(edge_params)
    a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  for step in 1:num_steps
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w3)
    if !isnothing(photon_params) && (step < num_steps)
      stochastic_radiation!(i, coords, photon_params..., ds_step)
    end
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  if !isnothing(edge_params)
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
  end
end


@makekernel fastgtpsa=true function order_eight_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, num_steps, edge_params, L)
  w0 =  1.7084530707869978*ds_step
  w1 =  0.102799849391985*ds_step
  w2 = -1.96061023297549*ds_step
  w3 =  1.93813913762276*ds_step
  w4 = -0.158240635368243*ds_step
  w5 = -1.44485223686048*ds_step
  w6 =  0.253693336566229*ds_step
  w7 =  0.914844246229740*ds_step
  if !isnothing(edge_params)
    a, tilde_m, Ksol, Kn0, e1, e2 = edge_params
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e1, 1)
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  for step in 1:num_steps
    ker(i, coords, params..., w7)
    ker(i, coords, params..., w6)
    ker(i, coords, params..., w5)
    ker(i, coords, params..., w4)
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w1)
    ker(i, coords, params..., w0)
    ker(i, coords, params..., w1) 
    ker(i, coords, params..., w2)
    ker(i, coords, params..., w3)
    ker(i, coords, params..., w4)
    ker(i, coords, params..., w5)
    ker(i, coords, params..., w6)
    ker(i, coords, params..., w7)
    if !isnothing(photon_params) && (step < num_steps)
      stochastic_radiation!(i, coords, photon_params..., ds_step)
    end
  end
  if !isnothing(photon_params)
    stochastic_radiation!(i, coords, photon_params..., ds_step / 2)
  end
  if !isnothing(edge_params)
    linear_bend_fringe!(i, coords, a, tilde_m, Ksol, Kn0, e2, -1)
  end
end