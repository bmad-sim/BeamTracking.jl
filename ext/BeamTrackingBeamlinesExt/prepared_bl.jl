# universal! builds the kernel chain of an element (slice) without launching it, as a
# PreKernelChain, and run! launches it. What depends on the Bunch at the call (bunch.t_ref) is kept
# in the chain as a PreParameter and evaluated by run! just before the launch. track!(bunch, bl)
# builds and runs each chain at every call; prepare builds them once, to track the same Bunch many
# times (e.g. turns of a ring) with track!(bunch, chains).

using BeamTracking: num_lower, time_lower, batch_lower
using StaticArrays: StaticArray

# g(bunch) is the value at the call. Like TimeDependentParam, it is a closure that calls the closures
# of the PreParameters it is computed from directly (a value computed from PreParameters is one):
# evaluating nested PreParameters through one recursive method would defeat type inference (a
# prepared ramp was 10x slower). Captured values must be concrete: a captured type is a DataType,
# and a captured variable assigned twice is boxed
struct PreParameter{G}
  g::G
end

# f(bunch, args...), where the arguments that are PreParameters are evaluated at the call
PreParameter(f, args::Tuple) = (let f = f, gs = map(_getter, args); PreParameter(bunch -> f(bunch, map(g -> g(bunch), gs)...)) end)
_getter(a::PreParameter) = a.g
_getter(a) = (let a = a; bunch -> a end)

process_preparameter(p::PreParameter, bunch) = p.g(bunch)
# Inside the arguments of a kernel (e.g. a tuple of parameters or a vector of strengths), only if
# the type contains a PreParameter
process_preparameter(a::Union{Tuple,NamedTuple,StaticArray}, bunch) =
  has_preparameter(typeof(a)) ? map(x -> process_preparameter(x, bunch), a) : a
process_preparameter(a, bunch) = a

# A compile-time constant (otherwise every argument of every chain is inspected at every call)
@generated has_preparameter(::Type{T}) where {T} = _has_preparameter(T)
_has_preparameter(T) = T <: PreParameter ||
  (T <: Union{Tuple,NamedTuple} && isconcretetype(T) && any(_has_preparameter, fieldtypes(T))) ||
  (T <: StaticArray && _has_preparameter(eltype(T)))

# Arithmetic as for TimeDependentParam: a value computed from a PreParameter is one, e.g. every
# strength computed from the reference momentum at the element when it is time dependent (without
# the per-particle ramp). Comparisons are errors: a decision that depends on the call cannot be
# made when preparing
for op in (:+, :-, :*, :/, :^)
  @eval begin
    Base.$op(a::PreParameter, b::PreParameter) = (let ga = a.g, gb = b.g; PreParameter(bunch -> $op(ga(bunch), gb(bunch))) end)
    Base.$op(a::PreParameter, b) = (let ga = a.g, b = b; PreParameter(bunch -> $op(ga(bunch), b)) end)
    Base.$op(a, b::PreParameter) = (let a = a, gb = b.g; PreParameter(bunch -> $op(a, gb(bunch))) end)
  end
end
# (Base has ^(x, ::Integer) for any x)
Base.:^(a::PreParameter, b::Integer) = (let ga = a.g, b = b; PreParameter(bunch -> ga(bunch)^b) end)
Base.literal_pow(::typeof(^), a::PreParameter, ::Val{N}) where {N} =
  (let ga = a.g; PreParameter(bunch -> Base.literal_pow(^, ga(bunch), Val{N}())) end)
for f in (:+, :-, :sqrt, :exp, :log, :sin, :cos, :tan, :sinh, :cosh, :tanh, :inv, :asin, :acos, :atan,
          :float, :abs, :sign, :conj)
  @eval Base.$f(a::PreParameter) = (let ga = a.g; PreParameter(bunch -> $f(ga(bunch))) end)
end
BeamTracking.atan2(a::PreParameter, b) = (let ga = a.g, b = b; PreParameter(bunch -> BeamTracking.atan2(ga(bunch), b)) end)
BeamTracking.atan2(a, b::PreParameter) = (let a = a, gb = b.g; PreParameter(bunch -> BeamTracking.atan2(a, gb(bunch))) end)
Base.broadcastable(a::PreParameter) = Ref(a)
for op in (:(<), :(<=), :isless, :(==), :isapprox)
  @eval Base.$op(::PreParameter, ::Any; kw...) = error("A value that depends on the call cannot be compared when preparing")
  @eval Base.$op(::Any, ::PreParameter; kw...) = error("A value that depends on the call cannot be compared when preparing")
end
# Lowered to the coordinates' number type at the call (as RefState fields and kernel arguments are)
BeamTracking.num_lower(::Type, a::PreParameter) = (let ga = a.g; PreParameter(b -> num_lower(eltype(b.coords.v), ga(b))) end)

# kc is the KernelChain, or a PreParameter whose value is the PreKernelChain (a chain built at the
# call, see universal!). After the launch, bunch.t_ref increases by dt, bunch.p_over_q_ref is
# p_over_q_ref_exit, and the position of the Bunch is s_exit and step_exit
struct PreKernelChain{K,D,P,S}
  kc::K
  dt::D
  p_over_q_ref_exit::P
  s_exit::S
  step_exit::Int
end

# The PreParameters are lowered as the RefState constructor, and make_kernel_call and push, lower
# what universal! passes there
# (also a tuple of parameters or a vector of strengths that contains PreParameters; the number type
# is that of the coordinates, which the RefState has, not captured: a captured type is a DataType)
process_field(a, bunch) = has_preparameter(typeof(a)) ? batch_lower(process_preparameter(a, bunch)) : a
process_arg(a, bunch) = has_preparameter(typeof(a)) ?
  num_lower(eltype(bunch.coords.v), time_lower(batch_lower(process_preparameter(a, bunch)))) : a

function process_preparameters(kc::KernelChain, bunch)
  T = BeamTracking.coordstype(kc.ref)
  call(c) = BeamTracking.KernelCall(c.kernel, map(a -> process_arg(a, bunch), c.args))
  r = kc.ref
  fields = (r.t_enter, r.beta_gamma_enter, r.t_exit, r.beta_gamma_exit, r.L, r.g, r.ds_step, r.s0, r.s_enter)
  ref = RefState{T}(map(a -> process_field(a, bunch), fields)..., r.step0, r.n_steps)
  return KernelChain(map(call, kc.chain), ref, kc.batch_start, map(call, kc.transforms_out),
                     map(call, kc.transforms_in))
end

# universal! launches its chain with this when tracking (one dynamic call per element)
function run_chain!(bunch::Bunch, kc::KernelChain; kwargs...)
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(bunch.coords, process_preparameters(kc, bunch); kwargs...)
  return nothing
end

# Not with run_chain!: a call with keyword arguments is not inlined, and p.kc passed to a call is
# copied out of p, which is boxed in the Vector of chains (each 8% slower for one particle on ESR)
function run!(bunch::Bunch, p::PreKernelChain; kwargs...)
  kc = process_preparameters(p.kc, bunch)
  dt = process_preparameter(p.dt, bunch)
  p_over_q_ref_exit = process_preparameter(p.p_over_q_ref_exit, bunch)
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(bunch.coords, kc; kwargs...)
  bunch.t_ref = bunch.t_ref + dt
  bunch.p_over_q_ref = p_over_q_ref_exit
  bunch.s = p.s_exit
  bunch.step = p.step_exit
  return nothing
end


# prepare: the chain is stored, and the next element is built for the reference momentum after
# this one, unless that is only known at the call (with a time dependent reference energy, where
# the next element only uses it at the call too)
# (with the position of bunch after the chain, which universal! has set)
function push_chain!(chains::Vector, bunch, kc, dt, p_over_q_ref_exit)
  push!(chains, PreKernelChain(kc, dt, p_over_q_ref_exit, bunch.s, bunch.step))
  p_over_q_ref_exit isa PreParameter || (bunch.p_over_q_ref = p_over_q_ref_exit)
  return nothing
end

# The chains of bl for bunch, which is not tracked. The keyword arguments are those of
# track!(bunch, bl) without the launch! ones
function prepare(
  bunch::Bunch,
  bl::Beamline;
  scalar_params::Bool=false,
  ramp_particle_energy_without_rf::Bool=false,
  ramp_update_each_particle::Bool=false,
  rf_on::Bool=true,
  batch_start=1,
  notify::Bool=true,
)
  chains = PreKernelChain[]
  length(bl.line) == 0 && return chains
  # Resolved once for the whole Beamline
  ref = resolve_reference!(bunch, Beamlines.trygetproperty(bl, :species_ref),
                           Beamlines.trygetproperty(bl, :p_over_q_ref), notify)
  context = bl.context
  # These are static for the whole Beamline, so the Vals are built once here and not per element
  v_ramp_particle_energy_without_rf = Val(ramp_particle_energy_without_rf)
  v_ramp_update_each_particle = Val(ramp_update_each_particle)
  v_shift_time = Val(bunch.species != ref.species_ref)
  # Built with a copy (at the beginning of the Beamline), whose reference momentum and position
  # push_chain! changes as the elements do
  b = Bunch(bunch.species, bunch.p_over_q_ref, bunch.t_ref, bunch.coords)
  for ele in bl.line
    @noinline _track!(b, ele, context, ref, scalar_params, v_ramp_particle_energy_without_rf,
                      v_ramp_update_each_particle, v_shift_time, rf_on, batch_start, chains)
  end
  return chains
end

function track!(bunch::Bunch, chains::Vector{PreKernelChain}; kwargs...)
  for pc in chains
    @noinline run!(bunch, pc; kwargs...)
  end
  return bunch
end
