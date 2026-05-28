@makekernel fastgtpsa=true function integrate_thin!(i, coords::Coords, mm, knl, ksl, a, g, tilde_m)
  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  end
  multipole_kick!(i, coords, mm, knl, ksl, -1)
  if !isnothing(coords.q)
    rotate_spin!(i, coords, a, g, tilde_m, mm, knl, ksl, 1/2)
  end
  execute_callbacks(coords, 0, (0, 0))
end