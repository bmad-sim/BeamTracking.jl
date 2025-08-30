@makekernel fastgtpsa=true function track_aperture_rectangular!(i, coords::Coords, x1, x2, y1, y2)
  v = coords.v

  coords.state[i] = vifelse(v[i,XI] < x1, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,XI] > x2, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,YI] < y1, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,YI] > y2, STATE_LOST, coords.state[i])
end

#

@makekernel fastgtpsa=true function track_aperture_elliptical!(i, coords::Coords, x1, x2, y1, y2)
  v = coords.v

  x = (2*v[i,XI] - (x2 + x1)) / (x2 - x1)
  y = (2*v[i,YI] - (y2 + y1)) / (y2 - y1)
  r = x*x + y*y
  coords.state[i] = vifelse(r > 1 || x2 < x1 || y2 < y1, STATE_LOST, coords.state[i])
end
