#=

Utility functions and "fake" APC. These will be moved to 
AcceleratorSimUtils.jl in the end.

=#



# Straight from SIMD.jl:
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.Vec{N, T}, v2::SIMD.ScalarTypes) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::SIMD.ScalarTypes, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::Bool, v1::T, v2::T) where {T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::SIMD.Vec{N, T}) where {N, T} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::T2, v2::SIMD.Vec{N, T}) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
@inline vifelse(v::SIMD.Vec{N, Bool}, v1::SIMD.Vec{N, T}, v2::T2) where {N, T, T2 <:SIMD.ScalarTypes} = SIMD.vifelse(v, v1, v2)
# Fallback for type unstable:
@inline vifelse(v::Union{Bool,SIMD.Vec{N, Bool}}, v1, v2) where {N} = ifelse(v, v1, v2)

#  Math =======================================================================
# u corresponds to unnormalized

# sinc/sincu is zero when the real part is Inf and imag is finite
isinf_real(x::Real) = isinf(x)
isinf_real(x::Number) = isinf(real(x)) && isfinite(imag(x))

# sinhc/sinhcu is zero when the imag part is Inf and real is finite
isinf_imag(x::Real) = false
isinf_imag(x::Number) = isfinite(real(x)) && isinf(imag(x))

# sincu copied from Boost library and correct limit behavior added
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinc.hpp
"""
    sincu(x)

Compute the unnormalized sinc function ``\\operatorname{sincu}(x) = \\sin(x) / (x)`` 
with accuracy near the origin.
"""
function sincu(x)
  #if isinf_real(x)
  #  return zero(x)
  #end

  threshold = 0.0004 # (120*eps(Float64))^(1/4)
  return vifelse(abs(x) > threshold, sin(x)/x, 1-x^2/6)
end

# sinhcu copied from Boost library and correct limit behavior added
# https://www.boost.org/doc/libs/1_87_1/boost/math/special_functions/sinhc.hpp

"""
    sinhcu(x)

Compute the unnormalized sinhc function ``\\operatorname{sinhcu}(x) = \\sinh(x) / (x)`` 
with accuracy near the origin.
"""
function sinhcu(x)
  #if isinf_imag(x)
  #  return zero(x)
  #end

  threshold = 0.0004 # (120*eps(Float64))^(1/4)
  return vifelse(abs(x) > threshold, sinh(x)/x, 1+x^2/6)
end


function atan2(y, x)
  return atan(y, x)
end


function atan2(y::SIMD.Vec{N, T}, x::SIMD.Vec{N, T}) where {N, T}
  arctan = atan(y/x)
  return vifelse(x > 0, arctan,
         vifelse((x < 0)  & (y >= 0),  arctan + SIMD.Vec{N, T}(T(pi)),
         vifelse((x < 0)  & (y < 0),   arctan - SIMD.Vec{N, T}(T(pi)),
         vifelse((x == 0) & (y > 0),   SIMD.Vec{N, T}(pi/2),
         vifelse((x == 0) & (y < 0),   SIMD.Vec{N, T}(-pi/2),
         vifelse((x == 0) & (y == 0),  SIMD.Vec{N, T}(0), 
         SIMD.Vec{N, T}(NaN)))))))
end


# Copy-pasted from sincc in bmad-ecosystem
function sincuc(x) 
  c0 = 1/6
  c1 = -1/120
  c2 = 1/5040
  c3 = -1/362880
  x2 = x^2
  return vifelse(abs(x) >= 0.1, (x-sin(x))/x^3, c0+x2*(c1+x2*(c2+x2*c3)))
end


"""
This function computes sin(sqrt(x))/sqrt(x) and cos(sqrt(x)), which are both 
necessary for exponentiating a rotation vector into a quaternion.
"""
function sincos_quaternion(x)
  threshold = 7.3e-8 # sqrt(24*eps(Float64))
  sq = sqrt(x)
  s, c = sincos(sq)
  s = s/sq
  s_out = vifelse(x > threshold, s, 1-x/6)
  c_out = vifelse(x > threshold, c, 1-x/2)
  return s_out, c_out
end


"""
This function computes sin(sqrt(x))/sqrt(x) and cos(sqrt(x)), which are both 
necessary for exponentiating a rotation vector into a quaternion.
"""
function sincos_quaternion(x::TPS{T}) where {T}
  ε = eps(T)
  N_max = 100
  N = 1
  conv_sin = false
  conv_cos = false
  y = one(x)
  prev_sin = one(x)
  prev_cos = one(x)
  result_sin = one(x)
  result_cos = one(x)
  #sq = one(x)
  # Using FastGTPSA! for the following makes other kernels run out of temps
  @FastGTPSA begin
    if x < 0.1
      while !(conv_sin && conv_cos) && N < N_max
        y = -y*x/((2*N)*(2*N - 1))
        result_sin = prev_sin + y/(2*N + 1)
        result_cos = prev_cos + y
        N += 1
        if normTPS(result_sin - prev_sin) < ε
          conv_sin = true
        end
        if normTPS(result_cos - prev_cos) < ε
          conv_cos = true
        end
        prev_sin = result_sin
        prev_cos = result_cos
      end
    else
      sq = sqrt(x)
      result_sin, result_cos = sincos(sq)
      result_sin = result_sin/sq
    end
  end
  if N == N_max
    @warn "sincos_quaternion convergence not reached in $N_max iterations"
  end
  return result_sin, result_cos
end


"""
This function computes exp(-i/2 v⋅σ) as a quaternion, where σ is the 
vector of Pauli matrices. If compute is false, it returns the identity quaternion.
"""
function expq(v, compute)
  n2 = @FastGTPSA (v[1]*v[1] + v[2]*v[2] + v[3]*v[3])/4
  n2_0 = zero(n2)
  s, c = sincos_quaternion(vifelse(compute, n2, n2_0))
  s = vifelse(compute, s, n2_0)
  return (c, s*v[1]/2, s*v[2]/2, s*v[3]/2)
end


#---------------------------------------------------------------------------------------------------
"""
This function computes exp(-i/2 v⋅σ) as a quaternion, where σ is the 
vector of Pauli matrices. If compute is false, it returns the identity quaternion.
"""
function quat_mul(q1, q20, q2x, q2y, q2z)
  """
  Returns q1 * q2.
  """
  a1, b1, c1, d1 = q1[Q0], q1[QX], q1[QY], q1[QZ]
  a2, b2, c2, d2 = q20, q2x, q2y, q2z
  @FastGTPSA begin
    a3 = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b3 = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c3 = a1*c2 - b1*d2 + c1*a2 + d1*b2
    d3 = a1*d2 + b1*c2 - c1*b2 + d1*a2
  end
  return (a3, b3, c3, d3)
end

#---------------------------------------------------------------------------------------------------

function quat_mul(q1, q2)
  return [
    q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3] - q1[4]*q2[4]
    q1[1]*q2[2] + q1[2]*q2[1] + q1[3]*q2[4] - q1[4]*q2[3]
    q1[1]*q2[3] - q1[2]*q2[4] + q1[3]*q2[1] + q1[4]*q2[2]
    q1[1]*q2[4] + q1[2]*q2[3] - q1[3]*q2[2] + q1[4]*q2[1]
  ]
end


#---------------------------------------------------------------------------------------------------

"""
    function quat_rotate!(r, q) -> rotated_r

Rotates vector `r` using quaternion `q`.
"""
@inline function quat_rotate!(r, q)

  w11 = 1 - 2*(q[QY]*q[QY] + q[QZ]*q[QZ])
  w12 =     2*(q[QX]*q[QY] - q[QZ]*q[Q0])
  w13 =     2*(q[QX]*q[QZ] + q[QY]*q[Q0])

  w21 =     2*(q[QX]*q[QY] + q[QZ]*q[Q0])
  w22 = 1 - 2*(q[QX]*q[QX] + q[QZ]*q[QZ])
  w23 =     2*(q[QY]*q[QZ] - q[QX]*q[Q0])

  w31 =     2*(q[QX]*q[QZ] - q[QY]*q[Q0])
  w32 =     2*(q[QY]*q[QZ] + q[QX]*q[Q0])
  w33 = 1 - 2*(q[QX]*q[QX] + q[QY]*q[QY])

  return (w11*r[1] + w12*r[2] + w13*r[3],
          w21*r[1] + w22*r[2] + w23*r[3],
          w31*r[1] + w32*r[2] + w33*r[3])
end


# Rotation matrix

"""
  rot_quaternion(x_rot, y_rot, z_rot)

Constructs a rotation quaternion based on the given Bryan-Tait angles.

Bmad/SciBmad follows the MAD convention of applying z, x, y rotations in that order.

The inverse quaternion reverses the order of operations and their signs.


Arguments:
- `x_rot::Number`: Rotation angle around the x-axis.
- `y_rot::Number`: Rotation angle around the y-axis.
- `z_rot::Number`: Rotation angle around the z-axis.

"""
function rot_quaternion(x_rot, y_rot, z_rot)
  qz = SA[cos(z_rot/2), 0, 0, sin(z_rot/2)]
  qx = SA[cos(x_rot/2), sin(x_rot/2), 0, 0]
  qy = SA[cos(y_rot/2), 0, sin(y_rot/2), 0]
  q = quat_mul(qx, qz[Q0], qz[QX], qz[QY], qz[QZ])
  q = quat_mul(qy, q[Q0], q[QX], q[QY], q[QZ])
  return SA[q[Q0], q[QX], q[QY], q[QZ]]
end

# Inverse rotation quaternion
function inv_rot_quaternion(x_rot, y_rot, z_rot)
  qz = SA[cos(z_rot/2), 0, 0, -sin(z_rot/2)]
  qx = SA[cos(x_rot/2), -sin(x_rot/2), 0, 0]
  qy = SA[cos(y_rot/2), 0, -sin(y_rot/2), 0]
  q = quat_mul(qx, qy[Q0], qy[QX], qy[QY], qy[QZ])
  q = quat_mul(qz, q[Q0], q[QX], q[QY], q[QZ])
  return SA[q[Q0], q[QX], q[QY], q[QZ]]
end

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
function compute_time(z, pz, ref)
  @FastGTPSA begin 
    K = (pz + 1)*ref.beta_gamma
    v = K/sqrt(1 + K*K)*C_LIGHT
    t = -z/v + ref.t
  end
  return t
end


#=
"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
@inline function bessel01_RF(x)
  threshold = 2.9e-7 # sqrt(64*eps(Float64))
  sq = sqrt(x)
  b0_out = besselj(0, sq)
  b1 = besselj(1, sq)
  b1 = b1/sq
  b1_out = ifelse(x > threshold, b1, 1/2-x/16)
  return b0_out, b1_out
end


"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
@inline function bessel01_RF(x::SIMD.Vec{N, T}) where {N, T}
  threshold = 2.9e-7 # sqrt(64*eps(Float64))
  sq = sqrt(x)
  b0_out = SIMDMathFunctions.vmap(besselj0, sq)
  b1 = SIMDMathFunctions.vmap(besselj1, sq)
  b1 = b1/sq
  b1_out = vifelse(x > threshold, b1, 1/2-x/16)
  return b0_out, b1_out
end


"""
This function computes J_0(sqrt(x)) and J_1(sqrt(x))/sqrt(x), which are 
necessary for tracking through a cylindrical pillbox cavity.
"""
function bessel01_RF(x::TPS{T}) where {T}
  ε = eps(T)
  N_max = 100
  N = 1
  conv0 = false
  conv1 = false
  y = one(x)
  prev0 = one(x)
  prev1 = one(x)
  result0 = one(x)
  result1 = one(x)/2
  @FastGTPSA begin
    while !(conv0 && conv1) && N < N_max
      y = -y*x/(4*N*N)
      result0 = result0 + y
      result1 = result1 + y/(2*(N + 1))
      N += 1
      if normTPS(result0 - prev0) < ε
        conv0 = true
      end
      if normTPS(result1 - prev1) < ε
        conv1 = true
      end
      prev0 = result0
      prev1 = result1
    end
  end
  if N == N_max
    @warn "bessel01_RF convergence not reached in $N_max iterations"
  end
  return result0, result1
end
=#

# Particle energy conversions =============================================================
R_to_E(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT*chargeof(species))^2 + massof(species)^2)
E_to_R(species::Species, E) = @FastGTPSA massof(species)*sinh(acosh(E/massof(species)))/C_LIGHT/chargeof(species) 
pc_to_R(species::Species, pc) = @FastGTPSA pc/C_LIGHT/chargeof(species)

R_to_gamma(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT/massof(species))^2+1)
R_to_pc(species::Species, R) = @FastGTPSA R*chargeof(species)*C_LIGHT
R_to_beta_gamma(species::Species, R) = @FastGTPSA R*chargeof(species)*C_LIGHT/massof(species)
R_to_v(species::Species, R) = @FastGTPSA chargeof(species)*C_LIGHT / sqrt(1+(massof(species)/(R*C_LIGHT))^2)
beta_gamma_to_v(beta_gamma) = @FastGTPSA C_LIGHT*beta_gamma/sqrt(1+beta_gamma^2)

#=


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
    R_ref(e_rest, beta_gamma, ne = 1)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the magnetic R_ref, ``B\\rho = p / q``.

DTA: Need to handle energy units other than ``\\mathrm{eV}``..
"""
function R_ref(e_rest, beta_gamma, ne = 1)
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
#function R_ref(e_rest, e_kin, ne = 1)
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

"""
    get_work(bunch::Bunch, ::Val{N}) where N -> work

Returns a tuple of `N` arrays of type `eltype(Bunch.v.x)` and 
length `length(Bunch.v.x)` which may be used as temporaries.

### Arguments
- `bunch`     -- Bunch to extract types and number of particles from
- `::Val{N}` -- Number of `N` temporary arrays desired
"""
function get_work(bunch::Bunch, ::Val{N}) where {N}
  sample = first(bunch.v.x)
  T = typeof(sample)
  N_particle = length(bunch.v.x)

  # Instead of using zeros, we do this to ensure 
  # same GTPSA descriptor if T isa TPS.
  return ntuple(Val{N}()) do t
    r = Vector{T}(undef, N_particle)
    for idx in eachindex(r)
      r[idx] = zero(sample)
    end
    r
  end
end

=#
