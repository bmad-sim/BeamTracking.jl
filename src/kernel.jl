# Get the register size for SIMD operations from VectorizationBase
const REGISTER_SIZE = VectorizationBase.register_size()

"""
    launch!(f!::F, v::V, args...; groupsize, multithread_threshold, use_KA, use_explicit_SIMD)

Launch a kernel function on particle coordinates with automatic optimization for both CPU and GPU backends.

# Arguments
- `f!`: Kernel function to execute. The kernel function `f!` must be of the form `f!(i, v, work, args...)`
- `v`: Input/output matrix of particle coordinates (always in SoA format)
- `args...`: Additional arguments for the kernel function

# Keyword Arguments
- `groupsize`: Number of threads per workgroup for GPU execution. If nothing, uses default based on register size for CPU
- `multithread_threshold`: Particle count threshold for enabling multithreading (default: 1750 * nthreads)
- `use_KA`: Whether to use KernelAbstractions.jl for execution (default: true for GPU, false for CPU with no groupsize)
- `use_explicit_SIMD`: Whether to use explicit SIMD vectorization (default: false)
"""
@inline function launch!(
  f!::F, 
  v::V, 
  args...; 
  groupsize::Union{Nothing,Integer}=nothing,
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=!(get_backend(v) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=false
) where {F<:Function,V}

  # Error handling
  # Cannot use both KA and explicit SIMD
  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end

  N_particle = size(v, 1)
  backend = get_backend(v)
  
  # GPU execution path
  if use_KA
    if !(backend isa GPU)
      error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
    end
    
    kernel! = isnothing(groupsize) ? f!(backend) : f!(backend, groupsize)
    kernel!(v, args...; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
    return v
  end

  # CPU execution path
  if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1
    execute_simd_cpu!(f!, v, N_particle, multithread_threshold, args...)
  else
    execute_standard_cpu!(f!, v, N_particle, multithread_threshold, args...)
  end
  
  return v
end

# Helper functions for CPU execution paths
@inline function execute_simd_cpu!(f!, v, N_particle, multithread_threshold, args...)
  # Get the SIMD lane width
  simd_lane_width = VectorizationBase.pick_vector_width(eltype(v))
  lane = VecRange{Int(simd_lane_width)}(0)
  # Calculate the number of SIMD-aligned particles
  rmn = rem(N_particle, simd_lane_width)
  N_SIMD = N_particle - rmn
  
  # Multithreaded SIMD
  if N_particle >= multithread_threshold
    Threads.@threads for i in 1:simd_lane_width:N_SIMD
      @assert last(i) <= N_particle "Out of bounds!"
      f!(lane+i, v, args...)
    end
  # Single-threaded SIMD
  else
    for i in 1:simd_lane_width:N_SIMD
      @assert last(i) <= N_particle "Out of bounds!"
      f!(lane+i, v, args...)
    end
  end
  
  # Process remaining particles
  for i in N_SIMD+1:N_particle
    @assert last(i) <= N_particle "Out of bounds!"
    f!(i, v, args...)
  end
end

@inline function execute_standard_cpu!(f!, v, N_particle, multithread_threshold, args...)
  # Multithreaded execution
  if N_particle >= multithread_threshold
    Threads.@threads for i in 1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, args...)
    end
  # Single-threaded execution with automatic vectorization
  else
    @simd for i in 1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, args...)
    end
  end
end

# TODO: collective effects
# May need to overload runkernel! for collective effects

"""
    runkernel!(f!::F, i, v, args...; kwargs...)

Execute a kernel either on a specific particle or a bunch of particles.

# Arguments
- `f!`: Kernel function to execute
- `i`: Particle index or nothing for a bunch
        If i is nothing, launches the kernel for a bunch with automatic optimization
        If i is an index, executes the kernel directly for that specific particle
- `v`: Input/output matrix of particle coordinates
- `args...`: Additional arguments for the kernel function
- `kwargs...`: Keyword arguments passed to launch! when executing in batch mode
"""
# When running kernel on a bunch, no index is provided, launch the kernel with automatic optimization
@inline runkernel!(f!::F, i::Nothing, v, args...; kwargs...) where {F} =launch!(f!, v, args...; kwargs...)
# When running kernel on a specific particle, execute the kernel directly for particle at that index
@inline runkernel!(f!::F, i, v, args...; kwargs...) where {F} = f!(i, v, args...)


"""
    @makekernel fcn

Macro to create a kernel function that can be executed on both CPU and GPU backends.
Transforms a regular function into a form compatible with KernelAbstractions.jl.

# Arguments
- `fcn`: Function definition to be transformed into a kernel

# Implementation Details
- Creates two versions of the function:
  1. A kernel version compatible with KernelAbstractions.jl
  2. The original function for direct CPU execution
- Handles const arguments appropriately for GPU execution
- Supports only positional arguments (no keyword arguments or default values)
"""
macro makekernel(fcn)
  fcn.head == :function || error("@makekernel must wrap a function definition")
  body = esc(fcn.args[2])
  signature = fcn.args[1].args

  fcn_name = esc(signature[1])
  args = esc.(signature[2:end])
  i = esc(signature[2])
  v = esc(signature[3])
  work = esc(signature[4])

  const_args = map(signature[5:end]) do t
    if t isa Expr
      if t.head == :(::)
        :(@Const($(esc(t))))
      else
        error("Default values and keyword arguments are NOT supported by @Const in KernelAbstractions.jl")
      end
    else
      :(@Const($(esc(t))))
    end
  end
  stripped_args = map(signature[2:end]) do t
    if t isa Expr
      if t.args[1] isa Expr
        esc(t.args[1].args[1])
      else
        esc(t.args[1])
      end
    else
      esc(t)
    end
  end

  return quote
    @kernel function $(fcn_name)($v, $work, $(const_args...))
      $(stripped_args[1]) = @index(Global, Linear)
      $(fcn_name)($(stripped_args...))
    end
  
    @inline function $(fcn_name)($(args...))
        $(body)
    end
  end
end