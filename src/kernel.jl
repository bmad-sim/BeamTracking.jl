
const REGISTER_SIZE = VectorizationBase.register_size()

# This is here in case kernel chain needs to be run 
# but is not fully filled. It does nothing
blank_kernel!(args...) = nothing

@kwdef struct KernelCall{K,A}
  kernel::K = blank_kernel!
  args::A   = ()
  function KernelCall(kernel, args)
    _args = map(t->time_lower(t), args)
    new{typeof(kernel),typeof(_args)}(kernel, _args)
  end 
end

# Store the state of the reference coordinate system
# Needed for time-dependent parameters
struct RefState{T,U}
  t::T          # Reference time
  beta_gamma::U # Reference energy
end

# Alias
struct KernelChain{C<:Tuple{Vararg{<:KernelCall}}, S<:Union{Nothing,RefState}}
  chain::C  # The tuple of KernelCalls
  ref::S    # An optional RefState for the initial time-dependent parameters
  KernelChain(chain, ref=nothing) = new{typeof(chain), typeof(ref)}(chain, ref)
end

KernelChain(::Val{N}, ref=nothing) where {N} = KernelChain(ntuple(t->KernelCall(), Val{N}()), ref)

push(kc::KernelChain, kcall::Nothing) = kc

push(kc::KernelChain, kcall) = @reset kc.chain = _push(kc.chain, kcall)

@unroll function _push(chain, kcall)
  i = 0
  @unroll for kcalli in chain
    i += 1
    if kcalli.kernel == blank_kernel!
      return @reset chain[i] = kcall
    end
  end
  error("Unable to push KernelCall to kernel chain: kernel chain is full")
end

# KA does not like Vararg
@kernel function generic_kernel!(coords::Coords, @Const(kc::KernelChain))
  i = @index(Global, Linear)
  @inline _generic_kernel!(i, coords, kc)
end

_generic_kernel!(i, coords, kc) = __generic_kernel!(i, coords, kc.chain, kc.ref)

@unroll function __generic_kernel!(i, coords::Coords, chain, ref)
  @unroll for kcall in chain
    args = process_args(i, coords, kcall.args, ref)
    (kcall.kernel)(i, coords, args...)
  end
  return nothing
end

function process_args(i, coords, args, ref)
  if !isnothing(ref) && static_timecheck(args) 
    let t = compute_time(coords.v[i,ZI], coords.v[i,PZI], ref)
      new_args = map(arg->teval(arg, t), args)
      return map(arg->teval(arg, t), args)
    end
  else
    return args
  end
end

# Generic function to launch a kernel on the bunch coordinates matrix
# Matrix v should ALWAYS be in SoA whether for real or as a view via tranpose(v)

@inline function launch!(
  coords::Coords{<:Any,V},
  kc::KernelChain;
  groupsize::Union{Nothing,Integer}=nothing, #backend isa CPU ? floor(Int,REGISTER_SIZE/sizeof(eltype(v))) : 256 
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=!(get_backend(coords.v) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=!use_KA # Default to use explicit SIMD on CPU
) where {V}
  v = coords.v
  N_particle = size(v, 1)

  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end
#=  
  if !use_KA && backend isa GPU
    error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
  end
=#
  if !use_KA
    if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1 # do SIMD
      simd_lane_width = VectorizationBase.pick_vector_width(eltype(V))
      lane = SIMD.VecRange{Int(simd_lane_width)}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because SIMD.VecRange SIMD
          _generic_kernel!(lane+i, coords, kc)
        end
      else
        for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because SIMD.VecRange SIMD
          _generic_kernel!(lane+i, coords, kc)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        _generic_kernel!(i, coords, kc)
      end
    else
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          _generic_kernel!(i, coords, kc)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          _generic_kernel!(i, coords, kc)
        end
      end
    end
  else
    backend = get_backend(v)
    if isnothing(groupsize)
      kernel! = generic_kernel!(backend)
    else
      kernel! = generic_kernel!(backend, groupsize)
    end
    kernel!(coords, kc; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
  end
  return v
end

# Call kernels directly
@inline runkernels!(i::Nothing, coords::Coords, kc::KernelChain; kwargs...) =  launch!(coords, kc; kwargs...)
@inline runkernels!(i, coords::Coords, kc::KernelChain; kwargs...) = _generic_kernel!(i, coords, kc)

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
@inline launch!(coords::Coords, kcall::KernelCall; kwargs...) = launch!(coords, KernelChain((kcall,)); kwargs...)

macro makekernel(args...)
  kwargs = args[1:length(args)-1]
  fcn = last(args)

  fcn.head == :function || error("@makekernel must wrap a function definition")
  body = esc(fcn.args[2])
  signature = fcn.args[1].args

  fcn_name = esc(signature[1])
  args = esc.(signature[2:end])

  # Check if function body contains a return:
  MacroTools.postwalk(body) do x
    !(@capture(x, return _)) || error("Return statement not permitted in a kernel function $(signature[1])")
  end

  check_kwargs(:makekernel, kwargs...)
  kwargnames = map(t->t[1], map(t->Pair(t.args...), kwargs))
  kwargvals = map(t->t[2],map(t->Pair(t.args...), kwargs))

  idx_fastgtpsa = findfirst(t->t==:fastgtpsa, kwargnames)
  idx_inbounds = findfirst(t->t==:inbounds, kwargnames)

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