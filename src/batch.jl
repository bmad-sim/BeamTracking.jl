#=

Given a batch = [k1, k2, k3], batch parameters are seen by the particles as

Particle 1: k1
Particle 2: k2
Particle 3: k3
Particle 4: k1
Particle 5: k2
Particle 6: k3

etc

=#

# BatchParam acts like a number. If V <: Number, then you 
# can operate arbitrary BatchParams together. If V <: AbstractArray, 
# then those arrays better be the same length, or else there
# will be an error.

# For BatchParams, we will have two types - one which is visible to the user 
# and can be operated on, and another which is done after "lowering" in a 
# KernelCall construction, which also stores the length of the batch array
# in the type itself. This is for constant folding + performance
struct BatchParam{V}
  batch::V
end

# BatchParam will act like a number
# Conversion of types to BatchParam
BatchParam(a::BatchParam) = a

# Make these apply via convert
Base.convert(::Type{D}, a::Number) where {D<:BatchParam} = BatchParam(a) # Scalar BatchParam
Base.convert(::Type{D}, a::BatchParam{T}) where {D<:BatchParam, T<:Number} = BatchParam(a) # Scalar BatchParam
Base.convert(::Type{BatchParam{T}}, a::D) where {T,D<:BatchParam} = BatchParam(convert(T, a.batch))

Base.zero(b::BatchParam) = BatchParam(zero(first(b.batch)))
Base.one(b::BatchParam)  = BatchParam(one(first(b.batch))) 

# Now define the math operations:
for op in (:+,:-,:*,:/,:^)
  @eval begin
    Base.$op(ba::BatchParam, b::Number)   = (let b = b; return BatchParam(map(x->$op(x, b), ba.batch)); end)
    Base.$op(a::Number,   bb::BatchParam) = (let a = a; return BatchParam(map(x->$op(a, x), bb.batch)); end)
    function Base.$op(ba::BatchParam, bb::BatchParam)
      if !(ba.batch isa Number) && !(bb.batch isa Number) && length(ba.batch) != length(bb.batch)
        error("Cannot perform operation $($op) with two non-scalar BatchParams of differing 
               lengths (received lengths $(length(ba.batch)) and $(length(ba.batch))).")
      end
      return BatchParam(map((x,y)->$op(x, y), ba.batch, bb.batch))
    end
  end
end

function Base.literal_pow(::typeof(^), ba::BatchParam, ::Val{N}) where {N}
  return BatchParam(map(x->Base.literal_pow(^, x, Val{N}()), ba.batch))
end

for t = (:+, :-, :sqrt, :exp, :log, :sin, :cos, :tan, :cot, :sinh, :cosh, :tanh, :inv,
  :coth, :asin, :acos, :atan, :acot, :asinh, :acosh, :atanh, :acoth, :sinc, :csc, :float,
  :csch, :acsc, :acsch, :sec, :sech, :asec, :asech, :conj, :log10, :isnan, :sign, :abs)
  @eval begin
    Base.$t(b::BatchParam) = BatchParam(map(x->($t)(x), b.batch))
  end
end

atan2(b1::BatchParam, b2::BatchParam) = BatchParam(map((x,y)->atan2(x,y), b1.batch, b2.batch))

for t = (:unit, :sincu, :sinhc, :sinhcu, :asinc, :asincu, :asinhc, :asinhcu, :erf, 
         :erfc, :erfcx, :erfi, :wf, :rect)
  @eval begin
    GTPSA.$t(b::BatchParam) = BatchParam(map(x->($t)(x), b.batch))
  end
end

Base.promote_rule(::Type{BatchParam{T}}, ::Type{U}) where {T, U<:Number}     = Base.promote_op(*, T, U)
Base.promote_rule(::Type{BatchParam{T}}, ::Type{BatchParam{S}}) where {T, S} = Base.promote_op(*, T, U)
Base.broadcastable(o::BatchParam) = Ref(o)

Base.isapprox(::TimeDependentParam, ::Number; kwargs...) = false
Base.isapprox(::Number, ::TimeDependentParam; kwargs...) = false
Base.:(==)(::TimeDependentParam, ::Number) = false
Base.:(==)(::Number, ::TimeDependentParam) = false
Base.isinf(::TimeDependentParam) = false


"""
    lane2vec(lane::SIMD.VecRange{N}) 
    
Given a SIMD.VecRange, will return an equivalent SIMD.Vec that
can be used in arithmetic operations for mapping integer indices
of particles to a given element in a batch.
"""
function lane2vec(lane::SIMD.VecRange{N}) where {N}
  # Try to match with vector register size, but 
  # only up to UInt32 -> ~4.3 billion particles, 
  # probably max on CPU...
  if Int(pick_vector_width(UInt32)) == N
    return SIMD.Vec{N,UInt32}(ntuple(i->lane.i+i-1, Val{N}()))
  else
    return SIMD.Vec{N,UInt64}(ntuple(i->lane.i+i-1, Val{N}()))
  end
end


struct _LoweredBatchParam{N,V}
  batch::V
  _LoweredBatchParam{N,V}(batch) where {N,V} = new{N,V}(batch)
end

function _LoweredBatchParam(bp::BatchParam)
  return _LoweredBatchParam{length(bp.batch),typeof(bp.batch)}(bp.batch)
end

# Necessary for GPU compatibility if batch is a GPU array
function Adapt.adapt_structure(to, lbp::_LoweredBatchParam{N}) where {N}
    batch = Adapt.adapt_structure(to, lbp.batch)
    return _LoweredBatchParam{N,typeof(batch)}(batch)
end

@inline teval(f::TimeFunction, t) = f(t)
@inline teval(f, t) = f
@inline teval(f::Tuple, t) = map(ti->teval(ti, t), f)

time_lower(tp::TimeDependentParam) = tp.f
time_lower(tp) = tp
time_lower(tp::Tuple) = map(ti->time_lower(ti), tp)
function time_lower(tp::SArray{N,TimeDependentParam}) where {N}
  f = Tuple(map(ti->ti.f, tp))
  return TimeFunction(t->SArray{N}(map(fi->fi(t), f)))
end
time_lower(tp::SArray{N,Any}) where {N} = time_lower(TimeDependentParam.(tp))

#static_timecheck(::Type{<:TimeFunction}) = true
static_timecheck(tp) = false
static_timecheck(::TimeFunction) = true
@unroll function static_timecheck(t::Tuple)
  @unroll for ti in t
    if static_timecheck(ti)
      return true
    end
  end
  return false
end
#static_timecheck(tp::T) where {T<:Tuple} = Val{any(t->static_timecheck(t) isa Val{true}, tp)}()