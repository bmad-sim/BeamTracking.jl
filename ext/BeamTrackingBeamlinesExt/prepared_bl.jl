# Prepared tracking: the kernel chain of every element (slice) of a Beamline is built once by
# universal!, which pushes it to a Vector (rec) instead of launching it, and is then launched at
# every call. What depends on the Bunch at the call (bunch.t_ref) is kept in the chain as a
# PreParameter and evaluated just before each launch.

using BeamTracking: num_lower, time_lower, batch_lower

# f(bunch, args...) is what universal! passes there at a call: a RefState field, a kernel argument,
# or what universal! sets after launching. f captures nothing (a captured variable can be boxed,
# and a captured type is stored as a DataType), so that it is type stable
struct PreParameter{F,A}
  f::F
  args::A
end

pre(::Nothing, f, bunch, args...) = f(bunch, args...)
pre(rec, f, bunch, args...) = PreParameter(f, args)

process_preparameter(p::PreParameter, bunch) = p.f(bunch, p.args...)
process_preparameter(a, bunch) = a

# The chain of a slice, and the bunch.t_ref and bunch.p_over_q_ref universal! sets after it
struct PreKernelChain{K,T,P}
  kc::K
  t_ref_exit::T
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

function _track!(bunch::Bunch, p::PreKernelChain; kwargs...)
  kc = process_preparameters(p.kc, bunch)
  t_ref_exit = process_preparameter(p.t_ref_exit, bunch)
  p_over_q_ref_exit = process_preparameter(p.p_over_q_ref_exit, bunch)
  @noinline launch!(bunch.coords, kc; kwargs...)
  bunch.t_ref = t_ref_exit
  bunch.p_over_q_ref = p_over_q_ref_exit
  return nothing
end

struct PreparedBeamline{V}
  chains::Vector{Any}
  species::Species
  p_over_q_ref::Any
end

# The keyword arguments are those of track!(bunch, bl) without the launch! ones. The Bunch is not
# tracked; a Bunch tracked with the result must have its species, p_over_q_ref and array type
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
  ref = resolve_reference!(bunch, Beamlines.trygetproperty(bl, :species_ref),
                           Beamlines.trygetproperty(bl, :p_over_q_ref), notify)
  context = bl.context
  v_ramp_particle_energy_without_rf = Val(ramp_particle_energy_without_rf)
  v_ramp_update_each_particle = Val(ramp_update_each_particle)
  v_shift_time = Val(bunch.species != ref.species_ref)
  # universal! sets the reference time and momentum of this copy, as when tracking
  b = Bunch(bunch.species, bunch.p_over_q_ref, bunch.t_ref, bunch.coords)
  rec = Any[]
  for ele in bl.line
    @noinline _track!(b, ele, context, ref, scalar_params, v_ramp_particle_energy_without_rf,
                      v_ramp_update_each_particle, v_shift_time, rf_on, batch_start, rec)
  end
  return PreparedBeamline{typeof(bunch.coords.v)}(rec, bunch.species, bunch.p_over_q_ref)
end

function track!(bunch::Bunch, pbl::PreparedBeamline{V}; kwargs...) where {V}
  same = bunch.species == pbl.species && isequal(bunch.p_over_q_ref, pbl.p_over_q_ref)
  same && bunch.coords.v isa V ||
    error("The Bunch differs from the one the Beamline was prepared with (species, " *
          "p_over_q_ref or coordinate array type)")
  for p in pbl.chains
    @noinline _track!(bunch, p; kwargs...)
  end
  return bunch
end
