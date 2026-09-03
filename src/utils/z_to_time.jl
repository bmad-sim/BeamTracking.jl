
# Particle z, pz to time
# Evaluate time-dependent arguments
# we need to get the particle time, for that we need particle's velocity
# We have pz = dP/P0 = (P-P0)/P0 = P/P0-1 = (gamma*beta)/(gamma0*beta0)-1
# so pz + 1 = (gamma*beta)/(gamma0*beta0) 
# And then
# (pz + 1)*beta_0*gamma_0 = gamma*beta = beta/sqrt(1-beta^2)
# [(pz + 1)*beta_0*gamma_0]^2*(1-beta^2) = beta^2
# [(pz + 1)*beta_0*gamma_0]^2 = beta^2*(1+[(pz + 1)*beta_0*gamma_0]^2)
# So
# beta = (pz + 1)*beta_0*gamma_0/sqrt(1+[(pz + 1)*beta_0*gamma_0]^2)
# 
# Therefore, we should pass to the kernel beta_0*gamma_0 and t_ref to get beta
function compute_time(z, pz, t_ref, beta_gamma_ref)
  @FastGTPSA begin 
    K = (pz + 1)*beta_gamma_ref
    v = K/sqrt(1 + K*K)*C_LIGHT
    t = -z/v + t_ref
  end
  return t
end