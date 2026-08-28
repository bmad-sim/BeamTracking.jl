
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
@generated function compute_time(z::T, pz, t_ref, beta_gamma_ref) where {T}
  if T == Float16 || T == Float32
    TC_LIGHT = T(C_LIGHT)
  else
    TC_LIGHT = C_LIGHT
  end
  return quote
    @FastGTPSA begin 
      K = (pz + 1)*beta_gamma_ref
      v = K/sqrt(1 + K*K)*$TC_LIGHT
      t = -z/v + t_ref
    end
    return t
  end
end