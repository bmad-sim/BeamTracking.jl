struct Map end 
# Map tracking method, not necessarily symplectic

@makekernel fastgtpsa=true function map!(i, coords::Coords, transport_map, transport_map_params, L)
  v = coords.v
  q = coords.q
  alive = (coords.state[i] == STATE_ALIVE)
  
  v_in = (v[i,XI], v[i,PXI], v[i,YI], v[i,PYI], v[i,ZI], v[i,PZI])

  if !isnothing(q)
    q_in = (q[i,Q0], q[i,QX], q[i,QY], q[i,QZ])
  else
    q_in = nothing
  end
  if isnothing(transport_map_params) # For non-breaking
    v_out, q_out = transport_map(v_in, q_in)
  else
    v_out, q_out = transport_map(v_in, q_in, transport_map_params)
  end
  v[i,XI]  = vifelse(alive, v_out[XI],  v[i,XI])
  v[i,PXI] = vifelse(alive, v_out[PXI], v[i,PXI])
  v[i,YI]  = vifelse(alive, v_out[YI],  v[i,YI])
  v[i,PYI] = vifelse(alive, v_out[PYI], v[i,PYI])
  v[i,ZI]  = vifelse(alive, v_out[ZI],  v[i,ZI])
  v[i,PZI] = vifelse(alive, v_out[PZI], v[i,PZI])

  if !isnothing(coords.q) && !isnothing(q_out)
    q[i,Q0] = vifelse(alive, q_out[Q0], q[i,Q0])
    q[i,QX] = vifelse(alive, q_out[QX], q[i,QX])
    q[i,QY] = vifelse(alive, q_out[QY], q[i,QY])
    q[i,QZ] = vifelse(alive, q_out[QZ], q[i,QZ])
  end
end

