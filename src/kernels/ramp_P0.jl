@makekernel fastgtpsa=true function update_P0!(i, coords::Coords, p_over_q_ref_initial, p_over_q_ref_final, ramp_without_rf)
  v = coords.v 
  v[i,PXI] = v[i,PXI] * p_over_q_ref_initial / p_over_q_ref_final
  v[i,PYI] = v[i,PYI] * p_over_q_ref_initial / p_over_q_ref_final
  if !ramp_without_rf
    v[i,PZI] = (p_over_q_ref_initial * (1 + v[i,PZI]) - p_over_q_ref_final) / p_over_q_ref_final
  end
end
