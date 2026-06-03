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
find_steps(::SciBmadStandard, L) = (1, L)
find_steps(::Any, L) = (1, L)