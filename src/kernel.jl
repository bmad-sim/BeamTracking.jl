
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
- `v`       -- Input/output matrix as an SoA or SoA view ALWAYS! (use transpose if AoS)
- `work`    -- A Vector of temporary vectors (columns of v) to run the kernel `f!`
- `args...` -- Any further arguments to run the kernel

## Keyword Arguments
- `simd_lane_width`       -- The number of SIMD lanes to use. Default is `REGISTER_SIZE/sizeof(eltype(A))`
- `multithread_threshold` -- Number of particles at which multithreading is used. Default is `1e6``
"""
@inline function launch!(
  f!::F, 
  idxs,
  args...; 
  simd::Val{simd_lane_width}=Val{0}(), # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=0 #Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {F<:Function,simd_lane_width}
  n = length(idxs)
  if idxs isa AbstractUnitRange && first(idxs) == 1 && simd_lane_width != 0 # do explicit SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(n, simd_lane_width)
    N_SIMD = n - rmn
    if n >= multithread_threshold
      Threads.@threads for i in range(start=1, stop=N_SIMD, step=simd_lane_width)
        f!(lane+i, args...)
      end
    else
      for i in range(start=1, stop=N_SIMD, step=simd_lane_width)
        f!(lane+i, args...)
      end
    end
    # Do the remainder
    for i in N_SIMD+1:n
      f!(i, args...)
    end
  else
    if n >= multithread_threshold
      Threads.@threads for i in idxs
        f!(i, args...)
      end
    else
      @simd for i in idxs
        f!(i, args...)
      end
    end
  end
  return nothing
end

# collective effects
# each threads corresponds to many particles
# go through each element, each thread loops through each 
# particle and does stuff with it

# Call launch!
@inline runkernel!(f!::F, idxs::AbstractArray, args...) where {F} = 
# JACC.parallel_for(last(idxs), f!, args...) 
#launch!(f!, idxs, args...)
#
#

# Call kernel directly
@inline runkernel!(f!::F, i::Integer, args...) where {F} = f!(i, args...)


#=

for particle in particles
  for ele in ring

  end
end

for ele in ring
  # do a bunch pre pro
  for particle in particle

  end
end
 =#