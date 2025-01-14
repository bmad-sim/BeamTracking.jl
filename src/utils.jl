#=

Utility functions needed for tracking.

=#


"""
    sr_gamma(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor ``\\gamma``.
"""
function sr_gamma(beta_gamma)
  return hypot(1, beta_gamma)
end



"""
    sr_gamma_m1(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor minus one, ``\\gamma - 1``.
"""
function sr_gamma_m1(beta_gamma)
  return beta_gamma^2 / (hypot(1, beta_gamma) + 1)
end


"""
    sr_beta(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the normalized velocity ``\\beta = v / c``.
"""
function sr_beta(beta_gamma)
  return beta_gamma / hypot(1, beta_gamma)
end


"""
    sr_pc(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the energy-equivalent momentum, ``pc``.
"""
function sr_pc(e_rest, beta_gamma)
  return e_rest * beta_gamma
end


"""
    sr_ekin(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the kinetic energy,
``E_\\text{kin} = mc^2(\\gamma - 1)``.
"""
function sr_ekin(e_rest, beta_gamma)
  return e_rest * sr_gamma_m1(beta_gamma)
end


"""
    sr_etot(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the total energy, ``E_\\tot = mc^2\\gamma``.
"""
function sr_etot(e_rest, beta_gamma)
  return e_rest * hypot(1, beta_gamma)
end


"""
    brho(e_rest, beta_gamma, ne = 1)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the magnetic rigidity, ``B\\rho = p / q``.

DTA: Need to handle energy units other than ``\\mathrm{eV}``..
"""
function brho(e_rest, beta_gamma, ne = 1)
  return (sr_pc(e_rest, beta_gamma) / (ne * C_LIGHT))
end
## If given ``E_\text{kin}`` instead of ``\beta\gamma``,
## use the following:
#
#function sr_gamma(e_rest, e_kin)
#  return e_kin / e_rest + 1
#end
#
#function sr_gamma_m1(e_rest, e_kin)
#  return e_kin / e_rest
#end
#
#function sr_beta_gamma(e_rest, e_kin)
#  return sqrt(e_kin / e_rest * (e_kin / e_rest + 2))
#end
#
#function sr_beta(e_rest, e_kin)
#  return sr_beta_gamma(e_rest, e_kin) / sr_gamma(e_rest, e_kin)
#end
#
#function sr_pc(e_rest, e_kin)
#  #return sqrt(e_kin * (e_kin + 2e_rest))
#  return e_rest * sr_beta_gamma(e_rest, e_kin)
#end
#
#function sr_etot(e_rest, e_kin)
#  return e_rest + e_kin
#end
#
#function brho(e_rest, e_kin, ne = 1)
#  return sr_pc(e_rest, e_kin) / (ne * clight)
#end

"""
    sincu(z)

## Description
Compute the unnormalized sinc function, ``\\operatorname{sincu}(z) = \\sin(z) / z``,
with a correct treatment of the removable singularity at the origin.

### Implementation
The function ``\\sin(z) / z = 1 - z^2 / 3! + z^4 / 5! - z^6 / 7! + \\cdots``.
For real values of ``z \\in (-1,1)``, one can truncate this series just before
the ``k^\\text{th}`` term, ``z^{2k} / (2k+1)!``, and the alternating nature of
the series guarantees the error does not exceed the value of this first truncated term.
It follows that if ``z^2 / 6 < \\varepsilon_\\text{machine}``, simply truncating
the series to 1 yields machine precision accuracy near the origin.
And outside that neighborhood, one simply computes ``\\sin(z) / z``.
On the otherhand, if we allow for complex values, one can no longer assume the
series alternates, and one must use a more conservative threshold.
Numerical exploration suggests that for ``|z| < 1`` the sum of terms starting
at the ``k^\\text{th}`` term is bounded by ``|z|^{2k} / (2k)!``.
It then follows that one should use ``|z|^2 / 2 < \\varepsilon_\\text{machine}``
as the criterion for deciding to truncate the series to 1 near the origin.
"""
function sincu(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sin(z) / z
  end
end

"""
    sinhcu(z)

## Description
Compute the hyperbolic sinc function, ``\\operatorname{sinhcu}(z) = \\operatorname{sinh}(z) / z``,
with a correct treatment of the removable singularity at the origin.

### Implementation
See sincu for notes about implementation.
"""
function sinhcu(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sinh(z) / z
  end
end

# These will go in GTPSA
sincu(z::TPS) = sinc(z/pi)
sinhcu(z::TPS) = sinhc(z/pi)
sincu(z::GTPSA.TempTPS) = (GTPSA.div!(z,z,pi); return GTPSA.__t_sinc(z);)
sinhcu(z::GTPSA.TempTPS) = (GTPSA.div!(z,z,pi); return GTPSA.__t_sinhc(z);)



