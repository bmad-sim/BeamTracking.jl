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


# === THIS BLOCK WAS PARTIALLY WRITTEN BY CLAUDE ===
# Generated function to unroll the ms, knl, ksl SArrays/tuples
# this is needed for type stability if the knl, ksl arrays are 
# NOT NTuples NOR SArrays, and have different eltypes. Also
# small performance benefit for unrolling and using for loops.
@generated function normalized_field(ms, knl, ksl, x, y, excluding)
    # Compute the promoted type at compile time
    N = length(ms) # ms will always be StaticArray - Int array of the indicies
    _process_arg(_t) = (_t <: SIMD.Vec ? eltype(_t) : _t)
    knltype = knl <: NTuple || knl <: SArray ? _process_arg(eltype(knl)) : promote_type(ntuple(i -> _process_arg(fieldtype(knl, i)), N)...)
    ksltype = ksl <: NTuple || ksl <: SArray ? _process_arg(eltype(ksl)) : promote_type(ntuple(i -> _process_arg(fieldtype(ksl, i)), N)...)
    T = promote_type(_process_arg(x), _process_arg(y), knltype, ksltype)
    quote
        knl_0 = zero($T)
        ksl_0 = zero($T)
        add = (ms[$N] != excluding && ms[$N] > 0)
        by = vifelse(add, knl[$N] * one($T), knl_0)
        bx = vifelse(add, ksl[$N] * one($T), ksl_0)

        $([quote
            for m in ms[$(j+1)]-1:-1:max(ms[$j], 1)
                t  = (by * x - bx * y) / m
                bx = (by * y + bx * x) / m
                by = t
                if m == ms[$j]
                  by += vifelse(ms[$j] != excluding, knl[$j] * one($T), knl_0)
                  bx += vifelse(ms[$j] != excluding, ksl[$j] * one($T), ksl_0)
                end
            end
        end for j in N-1:-1:1]...)

        return bx, by
    end
end
# === END CLAUDE ===



