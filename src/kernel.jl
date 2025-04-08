
const REGISTER_SIZE = VectorizationBase.register_size()
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6


# Generic function to launch a kernel on the bunch coordinates matrix
# Matrix v should ALWAYS be in SoA whether for real or as a view via tranpose(v)

"""
    launch!(f!::F, v, v0, work, args...; simd_lane_width, multithread_threshold)

General purpose function to launch a kernel `f!`. The syntax for a kernel `f!` must 
ALWAYS be the following:

## Arguments
- `i`       -- Particle index
- `v`       -- Output matrix as an SoA or SoA view ALWAYS! (use transpose if AoS)
- `v0`       -- Input matrix as an SoA or SoA view ALWAYS! (use transpose if AoS)
- `work`    -- A Vector of temporary vectors (columns of v) to run the kernel `f!`
- `args...` -- Any further arguments to run the kernel

## Keyword Arguments
- `simd_lane_width`       -- The number of SIMD lanes to use. Default is `REGISTER_SIZE/sizeof(eltype(A))`
- `multithread_threshold` -- Number of particles at which multithreading is used. Default is `1e6``
"""
@inline function launch!(
  f!::F, 
  v::A,
  v0::A,
  work, 
  args...; 
  simd_lane_width = floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {F,A}
  @assert size(v) == size(v0) "v and v0 must be the same size"
  N_particle = size(v, 1)
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, v0, work, args...)
      end
    else
      for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, v0, work, args...)
      end
    end
    # Do the remainder
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, v0, work, args...)
    end
  else
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, v0, work, args...)
      end
    else
      for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, v0, work, args...)
      end
    end
  end
  return v
end
