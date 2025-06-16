# Get the register size for SIMD operations from VectorizationBase
const REGISTER_SIZE = VectorizationBase.register_size()

# This is here in case kernel chain needs to be run 
# but is not fully filled. It does nothing
blank_kernel!(args...) = nothing

"""
    KernelCall{K,A}

A structure representing a single kernel call with its associated arguments.

# Fields
- `kernel::K`: The kernel function to be executed
- `args::A`: Tuple of arguments to be passed to the kernel
"""
@kwdef struct KernelCall{K,A}
  kernel::K = blank_kernel!
  args::A   = ()
end

# Alias - KernelChain is a tuple of KernelCalls
const KernelChain = Tuple{Vararg{<:KernelCall}}

"""
    push(kc::KernelChain, kcall)

Adds a new kernel call to the kernel chain at the first available empty slot.
Returns the modified kernel chain.

# Arguments
- `kc`: The kernel chain to modify
- `kcall`: The kernel call to add

# Returns
- Modified kernel chain with the new kernel call added

# Throws
- Error if the kernel chain is full
"""
@unroll function push(kc::KernelChain, kcall)
  i = 0
  @unroll for kcalli in kc
    i += 1
    # insert new kernel call at first empty slot
    if kcalli.kernel == blank_kernel!
      return @reset kc[i] = kcall
    end
  end
  error("Unable to push KernelCall to kernel chain: kernel chain is full")
end

# Constructors for KernelChain - create a chain of N blank kernel calls
KernelChain(::Val{N}) where {N} = ntuple(t->KernelCall(), Val{N}())
KernelChain(N::Integer) = ntuple(t->KernelCall(), Val{N}())

"""
    check_args(kc::KernelChain)

Validates the arguments of all kernel calls in the chain.
Currently does nothing.

# Arguments
- `kc`: The kernel chain to check

# Returns
- `true` if all arguments are valid
"""
@unroll function check_args(kc::KernelChain)
  @unroll for kcalli in kc
    check_args(kcalli)
  end
  return true
end

check_args(kcalli) = true

"""
    generic_kernel!(b::BunchView, kc::KernelChain)

The main kernel function that processes a bunch of particles.
Uses KernelAbstractions for parallel execution.

# Arguments
- `b`: The bunch view containing particle data
- `kc`: The kernel chain to execute
"""
@kernel function generic_kernel!(b::BunchView, @Const(kc::KernelChain))
  # Get the global index of the particle
  i = @index(Global, Linear)
  @inline _generic_kernel!(i, b, kc)
end

"""
    _generic_kernel!(i, b::BunchView, kc::KernelChain)

Internal function that executes all kernels in the chain for particle i.

# Arguments
- `i`: The particle index
- `b`: The bunchView containing particle data
- `kc`: The kernel chain to execute
"""
@unroll function _generic_kernel!(i, b::BunchView, kc::KernelChain)
  @unroll for kcall in kc
    (kcall.kernel)(i, b, kcall.args...)
  end
  return nothing
end

### Inline functions for launch! routing
# Launch kernel chain with explicit SIMD
@inline function _launch_simd(b::BunchView, kc::KernelChain, N_particle::Int, simd_lane_width::StaticInt, multithread_threshold::Int)
    lane = VecRange{Int(simd_lane_width)}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn

    if N_particle >= multithread_threshold
      # multi-thread
        Threads.@threads for i in 1:simd_lane_width:N_SIMD
            @assert last(i) <= N_particle "Out of bounds!"
            _generic_kernel!(lane+i, b, kc)
        end
    else
      # single-thread
        for i in 1:simd_lane_width:N_SIMD
            @assert last(i) <= N_particle "Out of bounds!"
            _generic_kernel!(lane+i, b, kc)
        end
    end
    # process particles not in SIMD chunks
    for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        _generic_kernel!(i, b, kc)
    end
end

# Launch kernel chain with standard CPU execution
@inline function _launch_cpu(b::BunchView, kc::KernelChain, N_particle::Int, multithread_threshold::Int)
    if N_particle >= multithread_threshold
      # multi-thread
        Threads.@threads for i in 1:N_particle
            @assert last(i) <= N_particle "Out of bounds!"
            _generic_kernel!(i, b, kc)
        end
    else
      # single-thread
        @simd for i in 1:N_particle
            @assert last(i) <= N_particle "Out of bounds!"
            _generic_kernel!(i, b, kc)
        end
    end
end

"""
    launch!(b::BunchView, kc::KernelChain; kwargs...)

Launches the kernel chain on the bunch.

# Arguments
- `b`: The bunch view containing particle data
- `kc`: The kernel chain to execute

# Keyword Arguments
- `groupsize`: Size of work groups for KernelAbstractions
- `multithread_threshold`: Particle count threshold for enabling multithreading
- `use_KA`: Whether to use KernelAbstractions
- `use_explicit_SIMD`: Whether to use explicit SIMD on CPU

# Returns
- The bunch after tracking
"""
@inline function launch!(
  b::BunchView{S,V},
  kc::KernelChain;
  groupsize::Union{Nothing,Integer}=nothing,
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=true,
  use_explicit_SIMD::Bool=!use_KA
) where {S,V}
  v = b.v
  N_particle = size(v, 1)
  backend = get_backend(v)

  # Validate configuration options
  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end
  if !use_KA && backend isa GPU
    error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
  end

  # Route for KernelAbstractions execution
  if use_KA
    if isnothing(groupsize)
      kernel! = generic_kernel!(backend)
    else
      kernel! = generic_kernel!(backend, groupsize)
    end
    # Launch kernel chain and wait for completion
    kernel!(b, kc; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
    return v
  end

  # Route for CPU execution without KA
  if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1
    # Launch kernel chain with explicit SIMD
    simd_lane_width = VectorizationBase.pick_vector_width(eltype(V))
    _launch_simd(b, kc, N_particle, simd_lane_width, multithread_threshold)
  else
    # Launch kernel chain with standard CPU execution
    _launch_cpu(b, kc, N_particle, multithread_threshold)
  end   
  return v
end

# Call kernels directly
@inline runkernels!(i::Nothing, b::BunchView, kc::KernelChain; kwargs...) =  launch!(b, kc; kwargs...)
@inline runkernels!(i, b::BunchView, kc::KernelChain; kwargs...) = _generic_kernel!(i, b, kc)

"""
    check_kwargs(mac, kwargs...)

Validates keyword arguments for kernel macros.

# Arguments
- `mac`: The macro name for error messages
- `kwargs`: The keyword arguments to validate

# Throws
- Error if invalid keyword arguments are provided
"""
function check_kwargs(mac, kwargs...)
  valid_kwargs = [:(fastgtpsa)=>Bool, :(inbounds)=>Bool]
  for k in kwargs
    if Meta.isexpr(k, :(=))
      pk = Pair(k.args...)
      idx = findfirst(t->t==pk[1], map(t->t[1], valid_kwargs))
      if isnothing(idx)
        error("Unrecognized input to @$(mac) macro: $(pk[1])")
      elseif typeof(pk[2]) != valid_kwargs[idx][2]
        error("Type for keyword argument `$(pk[1])` must be `$(valid_kwargs[idx][2])`")
      end
    else
      error("Unrecognized input to @$(mac) macro: $k")
    end
  end
end

# Also allow launch! on single KernelCalls
@inline launch!(b::BunchView, kcall::KernelCall; kwargs...) = launch!(b, (kcall,); kwargs...)

"""
    @makekernel

Macro for creating kernel functions with optional FastGTPSA and bounds checking.

# Optional keyword arguments:
  - `fastgtpsa`: Enable FastGTPSA optimization
  - `inbounds`: Enable bounds checking
"""
macro makekernel(args...)
  # Extract keyword arguments and function definition
  kwargs = args[1:length(args)-1]
  fcn = last(args)

  # Validate function definition format
  fcn.head == :function || error("@makekernel must wrap a function definition")
  body = esc(fcn.args[2])
  signature = fcn.args[1].args

  fcn_name = esc(signature[1])
  args = esc.(signature[2:end])

  # Check for return statements (not allowed in kernels)
  MacroTools.postwalk(body) do x
    !(@capture(x, return _)) || error("Return statement not permitted in a kernel function $(signature[1])")
  end

  # Validate keyword arguments
  check_kwargs(:makekernel, kwargs...)
  kwargnames = map(t->t[1], map(t->Pair(t.args...), kwargs))
  kwargvals = map(t->t[2],map(t->Pair(t.args...), kwargs))

  # Determine optimization flags
  idx_fastgtpsa = findfirst(t->t==:fastgtpsa, kwargnames)
  idx_inbounds = findfirst(t->t==:inbounds, kwargnames)

  # Generate appropriate function definition based on optimization flags
  if isnothing(idx_fastgtpsa) || !kwargvals[idx_fastgtpsa] # no fastgtpsa
    if isnothing(idx_inbounds) || kwargvals[idx_inbounds] # inbounds
      return quote
        @inline function $(fcn_name)($(args...))
          @inbounds begin
            $(body)
          end
        end
      end
    else # no inbounds
      return quote
        @inline function $(fcn_name)($(args...))
          $(body)
        end
      end
    end
  else # fastgtpsa
    if isnothing(idx_inbounds) || kwargvals[idx_inbounds] # inbounds
      return quote
        @inline function $(fcn_name)($(args...))
          @inbounds begin @FastGTPSA begin
            $(body)
          end end
        end
      end
    else # no inbounds
      return quote
        @inline function $(fcn_name)($(args...))
          @FastGTPSA begin
            $(body)
          end 
        end
      end
    end
  end
end