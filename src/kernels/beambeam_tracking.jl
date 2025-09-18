using ..GTPSA, ..BeamTracking, ..StaticArrays, ..ReferenceFrameRotations, ..KernelAbstractions, ..SIMD, ..SIMDMathFunctions, ..SpecialFunctions, ..Roots
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, STATE_ALIVE, STATE_LOST, @makekernel, Coords, BeamTracking.coord_rotation!
using ..BeamTracking: C_LIGHT

"""
	exact_beambeam!(i, coords::Coords, beta_0, gamsqr_0, tilde_m, 
					sig_x_strong, sig_y_strong, N_particles, n_slices,
					z_offset, beta_strong)

Track a particle through a beam-beam interaction with exact tracking.
The strong beam is modeled as n Gaussian slices.

## Arguments
- 'tilde_m'  -- mc2/p0c
- 'p0c'      -- reference momentum
- 'E_strong'    -- strong beam energy
- 'sig_x_strong':  horizontal std of strong beam
- 'sig_y_strong':  vertical std of strong beam
- 'sig_z_strong':  length of strong beam
- 'N_particles':   number of particles in strong bunch
- 'n_slices':      number of slices to divide strong beam into
- 'z_offset':      longitudinal offset of strong beam
"""
@makekernel function track_beambeam!(i, coords::Coords, p0c, E_strong, charge,
	sig_x_strong, sig_y_strong, sig_z_strong, N_particles,
	n_slices, z_offset)

	if(!isnothing(coords.q))
		error("Spin tracking not implemented with beam-beam yet")
	end

	v = coords.v

	rel_p = 1 + v[i, PZI]
	ps_02 = rel_p * rel_p - v[i, PXI] * v[i, PXI] - v[i, PYI] * v[i, PYI]
	good_momenta = (ps_02 > 0)
	alive_at_start = (coords.state[i] == STATE_ALIVE)
	coords.state[i] = vifelse(!good_momenta & alive_at_start, STATE_LOST, coords.state[i])
	alive = (coords.state[i] == STATE_ALIVE)

	pc = rel_p*p0c 
	#change to coord species later
	mc2 = BeamTracking.massof(Species("electron"))
	E_tot = sqrt(pc*pc + mc2*mc2)

	tilde_m = mc2 / p0c

	beta_0 = pc / E_tot
	gamsqr_0 = 1 / (1 - beta_0 * beta_0)

	mc2strong = BeamTracking.massof(Species("proton"))
	pc_strong = sqrt(E_strong*E_strong - mc2strong*mc2strong)

	beta_strong = pc_strong / E_strong

	r_e = 2.8179403262e-15
	bbi_const = -N_particles * charge * r_e / (2 * pi * p0c * (sig_x_strong + sig_y_strong))

	z_slices = bbi_slice_positions(n_slices, sig_z_strong) .- z_offset

	# For collision point calculation
	s0_factor = beta_0 / (beta_0 + beta_strong)

	# Begin at Ip
	s_lab = 0.0

	for slice_idx ∈ 1:n_slices
		z_slice = z_slices[slice_idx]

		slice_center = SA[0.0, 0.0, z_slice]

		s_lab_collision = find_collision_point(i, v, slice_center,
			beta_0, gamsqr_0, tilde_m, beta_strong,
			s_lab, s0_factor)

		del_s = s_lab_collision - s_lab

		if abs(del_s) > 0
			exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, del_s)
		end

		s_lab = s_lab_collision

		dx = v[i, XI] - slice_center[1]
		dy = v[i, YI] - slice_center[2]

		px_old = v[i, PXI]
		py_old = v[i, PYI]

		sigma = SA[sig_x_strong, sig_y_strong]
		nk_x, nk_y = bbi_kick_faddeeva(dx, dy, sigma)

		coef = bbi_const / n_slices
		kick_x = nk_x * coef
		kick_y = nk_y * coef

		new_px = px_old + kick_x
		new_py = py_old + kick_y
		v[i, PXI] = vifelse(alive, new_px, v[i, PXI])
		v[i, PYI] = vifelse(alive, new_py, v[i, PYI])

		e_factor = 0.25 / rel_p
		energy_change = e_factor * (kick_x * (kick_x + 2 * px_old) +
									kick_y * (kick_y + 2 * py_old))
		new_pz = v[i, PZI] + energy_change
		v[i, PZI] = vifelse(alive, new_pz, v[i, PZI])

		rel_p = 1 + v[i, PZI]
		new_beta = rel_p / sqrt(rel_p * rel_p + tilde_m * tilde_m)

	end
	exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, -s_lab)
end


"""
	find_collision_point(i, coords, v, slice_center, beta_0, gamsqr_0, tilde_m, 
						beta_strong, s_lab_current, s0_factor)

Find the s-position where the weak particle and strong beam slice meet using Brent's method.

## Arguments
- 'slice_center':  [x,y,z] position of strong beam slice
- 'beta_0':        Reference beta
- 'gamsqr_0':      Squared Lorentz factor
- 'tilde_m':       Normalized mass
- 'beta_strong':   Strong beam beta
- 's_lab_current': Current s-position in lab frame
- 's0_factor':     Velocity ratio factor

## Returns
- 's_lab':         Lab frame s-position of collision
"""
function find_collision_point(i, v, slice_center,
	beta_0, gamsqr_0, tilde_m, beta_strong,
	s_lab_current, s0_factor)
	@FastGTPSA begin
		# Save initial state
		v_save = SA[v[i, XI], v[i, PXI], v[i, YI], v[i, PYI], v[i, ZI], v[i, PZI]]

		function collision_func(s_lab_target)
			v[i, XI] = v_save[1]
			v[i, PXI] = v_save[2]
			v[i, YI] = v_save[3]
			v[i, PYI] = v_save[4]
			v[i, ZI] = v_save[5]
			v[i, PZI] = v_save[6]

			del_s = s_lab_target - s_lab_current

			if abs(del_s) > 0
				rel_p = 1 + v[i, PZI]
				P_t2 = v[i, PXI] * v[i, PXI] + v[i, PYI] * v[i, PYI]
				P_s2 = rel_p * rel_p - P_t2
				P_s = sqrt(max(P_s2, 1e-20))

				v[i, XI] += v[i, PXI] * del_s / P_s
				v[i, YI] += v[i, PYI] * del_s / P_s

				v[i, ZI] -= (rel_p * del_s * (P_t2 - v[i, PZI] * (2 + v[i, PZI]) / gamsqr_0) /
							 (beta_0 * sqrt(rel_p * rel_p + tilde_m * tilde_m) * P_s *
							  (beta_0 * sqrt(rel_p * rel_p + tilde_m * tilde_m) + P_s)))
			end

			s_particle = s_lab_target

			t_particle = s_lab_target / (beta_0 * C_LIGHT)
			s_slice = slice_center[3] + beta_strong * C_LIGHT * t_particle

			return s_particle - s_slice
		end

		# Initial guess
		s0 = slice_center[3] * s0_factor

		ds = max(abs(slice_center[3]) * 0.1, 0.001)

		try
			# Roots.find_zero with Brent method
			s_lab = find_zero(collision_func, (s0 - ds, s0 + ds), Roots.Brent(),
				atol = 1e-12, rtol = 1e-12)

			v[i, XI] = v_save[1]
			v[i, PXI] = v_save[2]
			v[i, YI] = v_save[3]
			v[i, PYI] = v_save[4]
			v[i, ZI] = v_save[5]
			v[i, PZI] = v_save[6]

			return s_lab
		catch
			# If fails, return initial guess
			v[i, XI] = v_save[1]
			v[i, PXI] = v_save[2]
			v[i, YI] = v_save[3]
			v[i, PYI] = v_save[4]
			v[i, ZI] = v_save[5]
			v[i, PZI] = v_save[6]
			return s0
		end
	end
end

"""
	bbi_slice_positions(n_slices, sig_z)

Calculate z-positions of beam-beam slices for a Gaussian distribution.

## Arguments
- 'n_slices': Number of slices
- 'sig_z':    Bunch length

## Returns
- Array of z-positions for each slice
"""
function bbi_slice_positions(n_slices, sig_z)
	if n_slices == 1
		return zeros(1) 
	end

	z_positions = zeros(n_slices) # fill with zeros
	for i ∈ 1:n_slices
		prob = (i - 0.5) / n_slices
		z_positions[i] = sig_z * sqrt(2) * erfinv(2 * prob - 1)
	end
	return z_positions
end

"""
	bbi_kick_faddeeva(dx, dy, sigma)

Beam-beam kick calculation using the Faddeeva complex error function.

## Arguments
- 'dx':      horizontal displacement from strong beam center
- 'dy':      vertical displacement from strong beam center
- 'sigma':   [sig_x, sig_y] Stds of strong beam

## Returns
- 'nk_x':    normalized horizontal kick
- 'nk_y':    normalized vertical kick
"""
function bbi_kick_faddeeva(dx, dy, sigma)
	sig_x, sig_y = sigma
	r = sig_y / sig_x
	x_norm = dx / sig_x
	y_norm = dy / sig_y

	# Round beam case
	if abs(r - 1.0) < 0.001

		amp = (x_norm*x_norm + y_norm*y_norm) / 2
		scale = 4 * pi * (1 - exp(-amp)) / amp
		scale = vifelse(amp > 30, 4 * pi / amp, scale)
		scale = vifelse(amp < 1e-4, 4 * pi, scale)

		return -x_norm * scale, -y_norm * scale
	end

	# Ensure r < 1 (swap if needed)
	flipped = (r > 1)
	if flipped
		r = 1 / r
		x_norm, y_norm = y_norm, x_norm
		sig_x, sig_y = sig_y, sig_x
	end

	denom = sqrt(2 * (1 - r^2))
	scale = 4 * sqrt(pi^3) * denom * (1 + r)
	u = abs(x_norm) / denom
	v = abs(y_norm) / denom
	sx = sign(x_norm)
	sy = sign(y_norm)

	z1 = complex(u, r * v)
	z2 = complex(r * u, v)

	# Faddeeva function: w(z) = exp(-z²) * erfc(-i*z)
	w1 = exp(-z1^2) * erfcx(-im * z1)
	w2 = exp(-z2^2) * erfcx(-im * z2)

	arg = (1 - r^2) * (u^2 + v^2)
	expon = exp(-arg)

	f = w1 - expon * w2

	nk_x = -scale * sx * imag(f)
	nk_y = -scale * sy * real(f)

	if flipped
		nk_x, nk_y = nk_y, nk_x
	end

	return nk_x, nk_y
end