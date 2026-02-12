@makekernel fastgtpsa=true function multipole_and_spin_kick!(i, coords, mm, knL, ksL, a, tilde_m, L)
  if isnothing(coords.q)
    multipole_kick!(i, coords, mm, knL, ksL, -1)
  else
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
    rotate_spin!(i, coords, a, 0, tilde_m, mm, knL./L, ksL./L, L)
    multipole_kick!(i, coords, mm, knL ./ 2, ksL ./ 2, -1)
  end
end
