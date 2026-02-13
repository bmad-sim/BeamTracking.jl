struct TimeFunction{F<:Function}
  f::F
end
(tf::TimeFunction)(t) = tf.f(t)

struct TimeDependentParam
  f::TimeFunction
  TimeDependentParam(f::TimeFunction=TimeFunction((t)->t)) = new(f)
  TimeDependentParam(f::Function) = new(TimeFunction(f))
end

# Convenience ctor
Time() = TimeDependentParam()

# Calling TimeDependentParam
(d::TimeDependentParam)(t) = d.f(t)

# Conversion of types to TimeDependentParam
TimeDependentParam(a::Number) = TimeDependentParam((t)->a) 
TimeDependentParam(a::TimeDependentParam) = a

# Make these apply via convert
Base.convert(::Type{D}, a::Number) where {D<:TimeDependentParam} = D(a)
Base.convert(::Type{D}, a::D) where {D<:TimeDependentParam} = a

Base.zero(::TimeDependentParam) = TimeDependentParam((t)->0)
Base.one(::TimeDependentParam) = TimeDependentParam((t)->1)

# Now define the math operations:
for op in (:+,:-,:*,:/,:^)
  @eval begin
    Base.$op(da::TimeDependentParam, b::Number)   = (let fa = da.f, b = b; return TimeDependentParam((t)-> $op(fa(t), b)); end)
    Base.$op(a::Number,   db::TimeDependentParam) = (let fb = db.f, a = a; return TimeDependentParam((t)-> $op(a, fb(t))); end)
    function Base.$op(da::TimeDependentParam, db::TimeDependentParam)
      let fa = da.f, fb = db.f
        return TimeDependentParam((t)-> $op(fa(t), fb(t)))
      end
    end
  end
end

function Base.literal_pow(::typeof(^), da::TimeDependentParam, ::Val{N}) where {N} 
  let fa = da.f
    return TimeDependentParam((t)->Base.literal_pow(^, fa(t), Val{N}()))
  end
end

for t = (:+, :-, :sqrt, :exp, :log, :sin, :cos, :tan, :cot, :sinh, :cosh, :tanh, :inv,
  :coth, :asin, :acos, :atan, :acot, :asinh, :acosh, :atanh, :acoth, :sinc, :csc, :float,
  :csch, :acsc, :acsch, :sec, :sech, :asec, :asech, :conj, :log10, :isnan, :sign, :abs)
  @eval begin
    Base.$t(d::TimeDependentParam) = (let f = d.f; return TimeDependentParam((t)-> ($t)(f(t))); end)
  end
end

atan2(d1::TimeDependentParam, d2::TimeDependentParam) = (let f1 = d1.f, f2 = d2.f; return TimeDependentParam((t)->atan2(f1(t),f2(t))); end)

for t = (:unit, :sincu, :sinhc, :sinhcu, :asinc, :asincu, :asinhc, :asinhcu, :erf, 
         :erfc, :erfcx, :erfi, :wf, :rect)
  @eval begin
    GTPSA.$t(d::TimeDependentParam) = (let f = d.f; return TimeDependentParam((t)-> ($t)(f(t))); end)
  end
end

Base.promote_rule(::Type{TimeDependentParam}, ::Type{U}) where {U<:Number} = TimeDependentParam
Base.broadcastable(o::TimeDependentParam) = Ref(o)

Base.isapprox(::TimeDependentParam, ::Number; kwargs...) = false
Base.isapprox(::Number, ::TimeDependentParam; kwargs...) = false
Base.:(==)(::TimeDependentParam, ::Number) = false
Base.:(==)(::Number, ::TimeDependentParam) = false
Base.isinf(::TimeDependentParam) = false

@inline teval(f::TimeFunction, t) = f(t)
@inline teval(f, t) = f

# === THIS BLOCK WAS PARTIALLY WRITTEN BY CLAUDE ===
# Generated function for arbitrary-length tuples
@generated function teval(f::T, t) where {T<:Tuple}
    N = length(T.parameters)
    # Use getfield with literal integer arguments
    exprs = [:(teval(Base.getfield(f, $i), t)) for i in 1:N]
    return :(tuple($(exprs...)))
end
# === END CLAUDE ===

time_lower(tp::TimeDependentParam) = tp.f
time_lower(tp) = tp
# We can use map on the CPU, but not the GPU. This step of time_lower-ing is on 
# the CPU and we are already type unstable here anyways, so we should do this.
time_lower(tp::T) where {T<:Tuple} = map(ti->time_lower(ti), tp)

# Arrays MUST be converted into tuples, for SIMD
time_lower(tp::SArray{N,TimeDependentParam}) where {N} = time_lower(Tuple(tp))

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