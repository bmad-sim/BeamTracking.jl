
@makekernel function multipole_and_spin_kick!(i, coords, mm, knL, ksL, a, tilde_m)
  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, knL, ksL, -1)
  else
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, knL, ksL, -1)
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
  end
end

@makekernel function multipole_and_spin_kick!(i, coords, mm, kn, ks, a, tilde_m, L)
  # (map: a batched length is a SIMD Vec, which can't be broadcast with a tuple)
  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, map(k -> k*L, kn), map(k -> k*L, ks), -1)
  else
    multipole_kick!(i, coords, mm, map(k -> k*(L/2), kn), map(k -> k*(L/2), ks), -1)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, kn, ks, -1, L)
    multipole_kick!(i, coords, mm, map(k -> k*(L/2), kn), map(k -> k*(L/2), ks), -1)
  end
end