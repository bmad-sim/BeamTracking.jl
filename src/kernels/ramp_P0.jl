"""
    reference_momentum_shift!(i, coords::Coords, P0c_old, dP0c, shift_pz)

Shift coordinates due to a change in reference energy `dE`.
"""
@inline function reference_momentum_shift!(i, coords::Coords, P0c_old, dP0c, ::Val{shift_pz}) where {shift_pz}
  @inbounds begin
    v = coords.v
    alive = (coords.state[i] == STATE_ALIVE)

    P0c_new = P0c_old + dP0c
    P0c_ratio = P0c_old / P0c_new
    
    v[i,PXI] = vifelse(alive, P0c_ratio * v[i,PXI], v[i,PXI])
    v[i,PYI] = vifelse(alive, P0c_ratio * v[i,PYI], v[i,PYI])
    if shift_pz
      v[i,PZI] = vifelse(alive, fma(P0c_old, v[i,PZI], -dP0c) / P0c_new, v[i,PZI])
    end
  end
end

# Note: Time() and TimeDependentParam always evaluate at the beginning of the element.
# e.g. Time dependent x_offset is the particle arrival time at the ENTRANCE of the element
# Time is NOT assumed to change in the element
# therefore, even the reference energy will be ramped uniformly to the bunch reference 
# time at the start of the element for ramp_update_each_particle
# In theory this could be changed in the future but is kind of a hassle and would be inconsistent
# with the rest of the framework and assumptions made for Time()

# However, to do this, we still need to know the initial time for the particles (for P0c_old)
# so we can't support it yet since it does require some work.
#=
@inline function callback_reference_momentum_shift!(i, coords, cur_s, cur_t_ref, P0c_old, dP0c, shift_pz)
  return reference_momentum_shift!(i, coords, P0c_old, dP0c, shift_pz)
end
=#

# Fallback
fma(a, b, c) = a * b + c