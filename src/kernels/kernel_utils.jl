
@makekernel function multipole_and_spin_kick!(i, coords, mm, knL, ksL, a, tilde_m)
  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, knL, ksL, -1)
  else
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, knL, ksL)
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
  end
end

@makekernel function multipole_and_spin_kick!(i, coords, mm, kn, ks, a, tilde_m, L)
  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, kn.*L, ks.*L, -1)
  else
    multipole_kick!(i, coords, mm, kn .* (L/2), ks .* (L/2), -1)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, kn, ks, L)
    multipole_kick!(i, coords, mm, kn .* (L/2), ks .* (L/2), -1)
  end
end
