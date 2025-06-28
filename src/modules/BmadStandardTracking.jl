struct BmadStandard end

module BmadStandardTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
using ..BeamTracking: C_LIGHT
const TRACKING_METHOD = BmadStandard

@makekernel fastgtpsa=true function linear_dipole_hard_edge_fringe!(i, b::BunchView, ge)
  v = b.v

  v[i,PXI] += ge * v[i,XI]                           # linear fringe
  v[i,PYI] -= ge * v[i,YI]
end

@makekernel fastgtpsa=true function combined_func!(i, b::BunchView, L, g, k0, k1, tilde_m)
  v = b.v

  rel_p  = 1 + v[i, PZI]                              # rel_p
  inv_rel_p = 1 / rel_p
  k_x  = k1 + g * k0                                  # k_x
  if k_x ≈ 0
    x_c = 0
  else
    x_c  = (g * rel_p - k0) / k_x                       # x_c
  end
  om_x  = sqrt(abs(k_x)) * sqrt(inv_rel_p)                      # om_x
  om_y  = sqrt(abs(k1)) * sqrt(inv_rel_p)                       # om_y

  arg = om_x * L
  if arg < 1e-6
    s_x = (1 - sign(k_x) * arg^2 / 6) * L      # s_x
    c_x = 1 - sign(k_x) * arg^2 / 2            # c_x
    z2 = g * L^2 / (2 * rel_p)                        # z2
  elseif k_x > 0
    s_x = sin(arg) / om_x                        # s_x
    c_x = cos(arg)                               # c_x
    z2 = -sign(k_x) * g * (1 - c_x) / (rel_p * om_x^2)# z2
  else
    s_x = sinh(arg) / om_x                       # s_x
    c_x = cosh(arg)                              # c_x
    z2 = -sign(k_x) * g * (1 - c_x) / (rel_p * om_x^2)# z2
  end

  arg = om_y * L
  if arg < 1e-6
    s_y = (1 + sign(k1) * arg^2 / 6) * L            # s_y
    c_y = 1 + sign(k1) * arg^2 / 2                  # c_y
  elseif k1 < 0
    s_y = sin(arg) / om_y                             # s_y
    c_y = cos(arg)                                    # c_y
  else
    s_y = sinh(arg) / om_y                            # s_y
    c_y = cosh(arg)                                   # c_y
  end

  # Save coordinates
  x_0  = v[i,  XI] - x_c                              # x = x - x_c
  px_0 = v[i, PXI]                                    # px
  y_0  = v[i,  YI]                                    # y
  py_0 = v[i, PYI]                                    # py

  # Update transverse
  v[i,  XI] = c_x * x_0 + s_x * px_0 * inv_rel_p + x_c
  v[i, PXI] = -sign(k_x) * om_x^2 * rel_p * s_x * x_0 + c_x * px_0
  v[i,  YI] = c_y * y_0 + s_y * py_0 * inv_rel_p
  v[i, PYI] = sign(k1) * om_y^2 * rel_p * s_y * y_0 + c_y * py_0

  # Longitudinal update
  # NEEDS CORRECTION
  # Update low_energy_z_correction for small pz case
  v[i, ZI] += L * (rel_p * sqrt(1 + tilde_m^2) / sqrt(rel_p^2 + tilde_m^2) - 1) + 
          (-g * x_c * L) +
          (-g * s_x) * x_0 +
          z2 * px_0 +
          (-sign(k_x) * om_x^2 * (L - c_x * s_x) / 4) * x_0^2 +
          (sign(k_x) * om_x^2 * s_x^2 / (2 * rel_p)) * x_0 * px_0 +
          (-(L + c_x * s_x) / (4 * rel_p^2)) * px_0^2 +
          (sign(k1) * om_y^2 * (L - c_y * s_y) / 4) * y_0^2 +
          (-sign(k1) * om_y^2 * s_y^2 / (2 * rel_p)) * y_0 * py_0 +
          (-(L + c_y * s_y) / (4 * rel_p^2)) * py_0^2
end

end