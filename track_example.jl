using Revise, BeamTracking, Beamlines, BenchmarkTools

function make_fodo(K1=0.40, L_quad=0.5, L_drift=5.0)
  qf = Quadrupole(K1=K1, L=L_quad, tracking_method=Linear())
  d1 = Drift(L=L_drift, tracking_method=Linear())
  qd = Quadrupole(K1=-qf.K1, L=L_quad, tracking_method=Linear())
  d2 = Drift(L=L_drift, tracking_method=Linear())
  return [qf, d1, qd, d2]
end

K1 = 0.40
L_quad = 0.5
L_drift = 5.0
N_fodo = 1000

bl = Beamline([ele for i in 1:N_fodo for ele in make_fodo(K1,L_quad,L_drift)]; Brho_ref=60.0)
N_particle = 10000

b0 = Bunch(N_particle, mem=BeamTracking.SoA)
@btime track!($b0, $bl)

bitbl = Beamlines.tobits(bl)
@btime track!($b0, $bitbl)