function find_steps(tm::BeamTracking.AbstractYoshida, L) 
  if L == 0
    return (1, L)
  end
  ds_step = tm.ds_step
  n_steps = tm.n_steps
  if ds_step < 0
    ds_step = L / n_steps
    return (n_steps, ds_step)
  else
    return (ceil(Int, L / ds_step), ds_step)
  end
end
find_steps(::Any, L) = (1, L)

# Temporary disgusting solution for callbacks - Yoshida
function compute_dt_ref(s, ker, params)
  beta_gamma_ref = find_beta_gamma_ref(ker, params)
  return s / beta_gamma_to_v(beta_gamma_ref)
end

find_beta_gamma_ref(::typeof(implicit_integrator!), params::P) where {P} = 1/params[3]
find_beta_gamma_ref(::typeof(exact_drift!), params::P) where {P} = 1/params[3]
find_beta_gamma_ref(::typeof(sks_multipole!), params::P) where {P} = 1/params[4]
find_beta_gamma_ref(::typeof(dkd_multipole!), params::P) where {P} = 1/params[4]
find_beta_gamma_ref(::typeof(bkb_multipole!), params::P) where {P} = 1/params[2]
find_beta_gamma_ref(::typeof(mkm_quadrupole!), params::P) where {P} = 1/params[4]
find_beta_gamma_ref(::typeof(exact_curved_drift!), params::P) where {P} = 1/params[7]
find_beta_gamma_ref(::K, params::P) where {K,P} = 0 # Fallback