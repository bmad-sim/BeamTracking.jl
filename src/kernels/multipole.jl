"""
    multipole_kick!(i, coords, ms, knl, ksl)

Track a beam of particles through a thin-lens multipole
having integrated normal and skew strengths listed in the
coefficient vectors knl and ksl respectively. The vector ms
lists the orders of the corresponding entries in knl and ksl.

The algorithm used in this function takes advantage of the
complex representation of the vector potential Az,
  - ``-Re{ sum_m (b_m + i a_m) (x + i y)^m / m! }``,
and uses a Horner-like scheme (see Shachinger and Talman
[SSC-52]) to compute the transverse kicks induced by a pure
multipole magnet. This method supposedly has good numerical
properties, though I've not seen a proof of that claim.

## Arguments
 - ms:  vector of m values for non-zero multipole coefficients
 - knl: vector of normal integrated multipole strengths
 - ksl: vector of skew integrated multipole strengths


     NB: Here the j-th component of knl (ksl) denotes the
       normal (skew) component of the multipole strength of
       order ms[j] (after scaling by the reference Bρ).
       For example, if ms[j] = 3, then knl[j] denotes the
       normal integrated sextupole strength scaled by Bρo.
       Moreover, and this is essential, the multipole
       coefficients must appear in ascending order.
"""
@makekernel fastgtpsa=true function multipole_kick!(i, coords::Coords, ms, knl, ksl, excluding)
  v = coords.v
  alive = (coords.state[i] == STATE_ALIVE)
  bx, by = normalized_field(ms, knl, ksl, v[i,XI], v[i,YI], excluding)
  bx_0 = zero(bx)
  by_0 = zero(by)
  v[i,PXI] -= vifelse(alive, by, by_0)                   
  v[i,PYI] += vifelse(alive, bx, bx_0)
end # function multipole_kick!()


function normalized_field(ms, knl, ksl, x, y, excluding)
  """
  Returns (bx, by), the transverse components of the magnetic field divided
  by the reference rigidty.
  """
  @FastGTPSA begin
    jm = length(ms)
    m  = ms[jm]
    add = (m != excluding && m > 0)
    knl_0 = zero(knl[jm]*x)
    ksl_0 = zero(ksl[jm]*y)
    by_0 = knl[jm]*one(x)
    bx_0 = ksl[jm]*one(y)
    by = vifelse(add, by_0, knl_0)
    bx = vifelse(add, bx_0, ksl_0)
    jm -= 1
    while 2 <= m
      m -= 1
      t  = (by * x - bx * y) / m
      bx = (by * y + bx * x) / m
      by = t
      add = (0 < jm && m == ms[jm])
      idx = max(1, jm) 
      new_by = by + vifelse(m != excluding, knl[idx], knl_0)
      new_bx = bx + vifelse(m != excluding, ksl[idx], ksl_0)
      by = vifelse(add, new_by, by)
      bx = vifelse(add, new_bx, bx)
      jm -= add
    end
  end
  return bx, by
end