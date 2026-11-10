# universal! builds the kernel chain of an element (slice) without launching it, as a
# PreKernelChain, and run! launches it. What depends on the Bunch at the call (bunch.t_ref) is kept
# in the chain as a PreParameter and evaluated by run! just before the launch. track!(bunch, bl)
# builds and runs each chain at every call; prepare builds them once, to track the same Bunch many
# times (e.g. turns of a ring) with track!(bunch, chains).

using BeamTracking: num_lower, time_lower, batch_lower

# f(bunch, args...) is the value at the call. f captures nothing (a captured variable can be boxed,
# and a captured type is stored as a DataType), so that it is type stable
struct PreParameter{F,A}
  f::F
  args::A
end

process_preparameter(p::PreParameter, bunch) = p.f(bunch, p.args...)
process_preparameter(a, bunch) = a

# kc is the KernelChain, or a PreParameter whose value is the PreKernelChain (a chain built at the
# call, see universal!). After the launch, bunch.t_ref increases by dt and bunch.p_over_q_ref is
# p_over_q_ref_exit
struct PreKernelChain{K,D,P}
  kc::K
  dt::D
  p_over_q_ref_exit::P
end

# The PreParameters are lowered as the RefState constructor, and make_kernel_call and push, lower
# what universal! passes there
function process_preparameters(kc::KernelChain, bunch)
  T = BeamTracking.coordstype(kc.ref)
  field(a) = a isa PreParameter ? batch_lower(process_preparameter(a, bunch)) : a
  arg(a) = a isa PreParameter ?
    num_lower(T, time_lower(batch_lower(process_preparameter(a, bunch)))) : a
  call(c) = BeamTracking.KernelCall(c.kernel, map(arg, c.args))
  r = kc.ref
  fields = (r.t_enter, r.beta_gamma_enter, r.t_exit, r.beta_gamma_exit, r.L, r.g, r.ds_step, r.s0)
  ref = RefState{T}(map(field, fields)...)
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
  # noinline necessary here for small binaries and faster execution
  @noinline launch!(bunch.coords, kc; kwargs...)
  bunch.t_ref = bunch.t_ref + p.dt
  bunch.p_over_q_ref = p.p_over_q_ref_exit
  return nothing
end

run!(bunch::Bunch, p::PreKernelChain{<:PreParameter}; kwargs...) =
  run!(bunch, process_preparameter(p.kc, bunch); kwargs...)

# track!: a chain built at the call is run at once (the others are run by universal! itself)
push_chain!(target::Bunch, bunch, pc; kwargs...) = run!(target, pc; kwargs...)

# A chain built at the call (see universal!) is returned to run!
push_chain!(::Nothing, bunch, pc; kwargs...) = pc

# prepare: the chain is stored, and the next element is built for the reference momentum after
# this one, unless this one is built at the call (which leaves the reference momentum unchanged
# when it is not time dependent)
function push_chain!(chains::Vector, bunch, pc)
  push!(chains, pc)
  isnothing(pc.p_over_q_ref_exit) || (bunch.p_over_q_ref = pc.p_over_q_ref_exit)
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
  # Built with a copy, whose reference momentum push_chain! changes as the elements do
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
