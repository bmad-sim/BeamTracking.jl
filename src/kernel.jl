
const REGISTER_SIZE = VectorizationBase.register_size()
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6

@kwdef struct KernelCall{K,A}
  kernel::K = blank_kernel!
  args::A   = ()
end

# In theory one can chain an entire lattice together - but that creates a giant type (N elements 
# type parameters) and this puts enormous strain on the compiler, to the point where it is actually
# slower (I tested it using tuples/@generated functions before )

# So this KernelChain allows a balance between giant method lookup table + stress on the compiler, 
# and some optimizations by chaining KernelCalls together
@kwdef struct KernelChain{K1<:KernelCall,K2<:KernelCall,K3<:KernelCall,K4<:KernelCall,K5<:KernelCall}
  k1::K1 = KernelCall()
  k2::K2 = KernelCall()
  k3::K3 = KernelCall()
  k4::K4 = KernelCall()
  k5::K5 = KernelCall()
end

function push(kc, kcall)
  if kc.k1.kernel == blank_kernel!
    return @reset kc.k1 = kcall
  elseif kc.k2.kernel == blank_kernel!
    return @reset kc.k2 = kcall
  elseif kc.k3.kernel == blank_kernel!
    return @reset kc.k3 = kcall
  elseif kc.k4.kernel == blank_kernel!
    return @reset kc.k4 = kcall
  elseif kc.k5.kernel == blank_kernel!
    return @reset kc.k5 = kcall
  else
    error("Maximum allowed length of KernelChain (5) exceeded!")
  end
end

# KA does not like Vararg
@kernel function generic_kernel!(@Const(kc::KernelChain), com_args)
  i = @index(Global, Linear)
  @inline generic_kernel!(i, kc, com_args...)
end

@inline function generic_kernel!(i, kc::KernelChain, com_args...)
  (kc.k1.kernel)(i, com_args..., kc.k1.args...)
  (kc.k2.kernel)(i, com_args..., kc.k2.args...)
  (kc.k3.kernel)(i, com_args..., kc.k3.args...)
  (kc.k4.kernel)(i, com_args..., kc.k4.args...)
  (kc.k5.kernel)(i, com_args..., kc.k5.args...)
  return nothing
end

blank_kernel!(args...) = nothing

@inline function launch!(
  kc::KernelChain,
  com_args::Vararg{V};
  groupsize::Union{Nothing,Integer}=nothing, #backend isa CPU ? floor(Int,REGISTER_SIZE/sizeof(eltype(v))) : 256 
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=!(get_backend(first(com_args)) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=false
) where {V}
  v = first(com_args)
  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end
  N_particle = size(v, 1)
  backend = get_backend(v)
  if !use_KA && backend isa GPU
    error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
  end

  if !use_KA
    if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1 # do SIMD
      simd_lane_width = VectorizationBase.pick_vector_width(eltype(V))
      lane = VecRange{Int(simd_lane_width)}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          generic_kernel!(lane+i, kc, com_args...)
        end
      else
        for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          generic_kernel!(lane+i, kc, com_args...)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        generic_kernel!(i, kc, com_args...)
      end
    else
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          generic_kernel!(i, kc, com_args...)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          generic_kernel!(i, kc, com_args...)
        end
      end
    end
  else
    if isnothing(groupsize)
      kernel! = generic_kernel!(backend)
    else
      kernel! = generic_kernel!(backend, groupsize)
    end
    kernel!(kc, com_args; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
  end
  return v
end



# Generic function to launch a kernel on the bunch coordinates matrix
# Matrix v should ALWAYS be in SoA whether for real or as a view via tranpose(v)

#=
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
  v::V, 
  args...; 
  groupsize::Union{Nothing,Integer}=nothing, #backend isa CPU ? floor(Int,REGISTER_SIZE/sizeof(eltype(v))) : 256 
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=!(get_backend(v) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=false
) where {F<:Function,V}

  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end

  N_particle = size(v, 1)
  backend = get_backend(v)
  if !use_KA && backend isa GPU
    error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
  end

  if !use_KA
    if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1 # do SIMD
      simd_lane_width = VectorizationBase.pick_vector_width(eltype(V))
      lane = VecRange{Int(simd_lane_width)}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          f!(lane+i, v, args...)
        end
      else
        for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          f!(lane+i, v, args...)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, args...)
      end
    else
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          f!(i, v, args...)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          f!(i, v, args...)
        end
      end
    end
  else
    if isnothing(groupsize)
      kernel! = f!(backend)
    else
      kernel! = f!(backend, groupsize)
    end
    kernel!(v, args...; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
  end
  return v
end
=#
# collective effects
# each threads corresponds to many particles
# go through each element, each thread loops through each 
# particle and does stuff with it

# Call launch!
#@inline runkernel!(f!::F, i::Nothing, v, args...; kwargs...) where {F} =launch!(f!, v, args...; kwargs...)
@inline runkernels!(i::Nothing, kc::KernelChain, com_args...; kwargs...) =  launch!(kc, com_args...; kwargs...)

# Call kernels directly
@inline runkernels!(i, kc::KernelChain, com_args...; kwargs...) = generic_kernel!(i, kc, com_args...)

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

macro makekernel(args...)
  kwargs = args[1:length(args)-1]
  fcn = last(args)

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

  
  check_kwargs(:makekernel, kwargs...)
  kwargnames = map(t->t[1], map(t->Pair(t.args...), kwargs))
  kwargvals = map(t->t[2],map(t->Pair(t.args...), kwargs))

  idx_fastgtpsa = findfirst(t->t==:fastgtpsa, kwargnames)
  idx_inbounds = findfirst(t->t==:inbounds, kwargnames)

  if isnothing(idx_fastgtpsa) || !kwargvals[idx_fastgtpsa] # no fastgtpsa
    if isnothing(idx_inbounds) || kwargvals[idx_inbounds] # inbounds
      return quote
        @kernel function $(fcn_name)($v, $work, $(const_args...))
          $(stripped_args[1]) = @index(Global, Linear)
          $(fcn_name)($(stripped_args...))
        end
      
        @inline function $(fcn_name)($(args...))
          @inbounds begin
            $(body)
          end
        end
      end
    else # no inbounds
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
  else # fastgtpsa
    if isnothing(idx_inbounds) || kwargvals[idx_inbounds] # inbounds
      return quote
        @kernel function $(fcn_name)($v, $work, $(const_args...))
          $(stripped_args[1]) = @index(Global, Linear)
          $(fcn_name)($(stripped_args...))
        end
      
        @inline function $(fcn_name)($(args...))
          @inbounds begin @FastGTPSA! begin
            $(body)
          end end
        end
      end
    else # no inbounds
      return quote
        @kernel function $(fcn_name)($v, $work, $(const_args...))
          $(stripped_args[1]) = @index(Global, Linear)
          $(fcn_name)($(stripped_args...))
        end
      
        @inline function $(fcn_name)($(args...))
          @FastGTPSA! begin
            $(body)
          end 
        end
      end

    end
  end

end


macro localvars(work, vars_and_block...)
  vars = vars_and_block[1:end-1]
  block = last(vars_and_block)
  all(t->t isa Symbol, vars) || error("Invalid input for @localvars")
  block = MacroTools.postwalk(esc(block)) do x
    for i in 1:length(vars)
      if x == vars[i]
        worki = MacroTools.postwalk(work) do x
          if x == :_
            return :($i)
          else
            return x
          end
        end
        return :($(worki))
      end
    end
    return x
  end
end

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