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
- 'beta_a_strong' -- strong beam beta in plane a 
- 'alpha_a_strong' -- strong beam alpha in plane a
- 'beta_b_strong' -- strong beam beta in plane b 
- 'alpha_b_strong' -- strong beam alpha in plane b
"""
@makekernel function track_beambeam!(i, coords::Coords, p0c, E_strong, charge,
	sig_x_strong, sig_y_strong, sig_z_strong, N_particles,
	n_slices,x_offset,y_offset,x_pitch,y_pitch, z_offset,tilt, 
	beta_a_strong, alpha_a_strong,
	beta_b_strong, alpha_b_strong, crab, crab_tilt)

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
    print("tilde_m: ", tilde_m, "\n")

	beta_0 = pc / E_tot
	gamsqr_0 = 1 / (1 - beta_0 * beta_0)

	part_time1 = v[i, ZI] / (beta_0 * C_LIGHT)
	
	mc2strong = BeamTracking.massof(Species("proton"))
	pc_strong = sqrt(E_strong*E_strong - mc2strong*mc2strong)

	beta_strong = pc_strong / E_strong

	r_e = 1.4399645468825422e-9
	bbi_const = -N_particles * charge * r_e / (2 * pi * p0c * (sig_x_strong + sig_y_strong))

	z_slices = bbi_slice_positions(n_slices, sig_z_strong)

	# For collision point calculation
    #check if tps TODO
	s0_factor = beta_0 / (beta_0 + beta_strong)

	# Begin at Ip
    part_time0 = -v[i,ZI]/(beta_0*C_LIGHT)
    if length(v[i,ZI]) > 1
	    s_lab = SIMD.Vec(zeros(length(v[i,ZI]))...)
        time = Ref(part_time0)
    else
        s_lab = 0.0
        time = Ref(part_time0)
    end
	# v[i,:] = offset_particle(i, v, x_offset, y_offset, 
    #                     z_offset, x_pitch, y_pitch, 
    #                     tilt, true)

	sigmaini = SA[sig_x_strong, sig_y_strong]
    s00 = (z_offset + v[i,ZI]) * s0_factor
    println("v[i,XI] before beam-beam: ", v[i,XI], "\n")
	for slice_idx ∈ 1:n_slices
		z_slice = z_slices[slice_idx]
        # print("coords before drift: X - ", coords.v[i,XI],", PX - ", coords.v[i,PXI],
        #       ", Y - ", coords.v[i,YI],", PY - ", coords.v[i,PYI],
        #       ", Z - ", coords.v[i,ZI],", PZ - ", coords.v[i,PZI],"\n")
        # exact_multi_drift!(i, coords, beta_0, gamsqr_0, tilde_m, 10.0)
        # print("coords after drift: X - ", coords.v[i,XI],", PX - ", coords.v[i,PXI],
        #       ", Y - ", coords.v[i,YI],", PY - ", coords.v[i,PYI],
        #       ", Z - ", coords.v[i,ZI],", PZ - ", coords.v[i,PZI],"\n")

        # exact_multi_drift!(i, coords, beta_0, gamsqr_0, tilde_m, -10)

		slice_center = strong_beam_center(z_slice, crab, crab_tilt)
		s_lab_collision = find_collision_point(i, coords, v, slice_center, 
                                              beta_0, gamsqr_0, tilde_m, beta_strong,
                                              s_lab, s0_factor,s00,z_offset, part_time1, p0c,time,part_time0)
        println("v[i,XI] before collision: ", v[i,XI], "\n")
        print("s_lab_collision: ", s_lab_collision, "\n")
        del_s = s_lab_collision - s_lab
        p_rel = 1 + v[i,PZI]
        px_rel = v[i,PXI]/p_rel
        py_rel = v[i,PYI]/p_rel
        pxy2 = px_rel*px_rel + py_rel*py_rel
        ps_rel = sqrt(1 - pxy2)
        dt = del_s / (beta_0 * C_LIGHT* ps_rel)
        time[] = time[] + dt
        print("time after col: ", time[], "\n")
		exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, del_s)
        println("v[i,XI] after drift to collision: ", v[i,XI], "\n")

		s_lab = s_lab_collision

		dx = v[i, XI] - slice_center[1]
		dy = v[i, YI] - slice_center[2]

		px_old = v[i, PXI]
		py_old = v[i, PYI]
		sigma, bbi_const, dsigma_ds = strong_beam_sigma_calc(sigmaini[1], sigmaini[2], 
                                s_lab,
                                beta_a_strong, alpha_a_strong,
                                beta_b_strong, alpha_b_strong,
                                N_particles * charge * r_e, 
                                p0c)
		nk_x, nk_y, dnk_unscaled = bbi_kick_faddeeva(dx, dy, sigma)
		coef = bbi_const / n_slices
        print("coef: ", coef, "\n")
		kick_x = nk_x * coef
        print("kick_x: ", kick_x, "\n")
		kick_y = nk_y * coef
        print("kick_y: ", kick_y, "\n")
		dnk = (
                (coef * dnk_unscaled[1][1], coef * dnk_unscaled[1][2]),
                (coef * dnk_unscaled[2][1], coef * dnk_unscaled[2][2])
            )
        print("dnk: ", dnk, "\n")

		new_px = px_old + kick_x
		new_py = py_old + kick_y

		v[i, PXI] = vifelse(alive, new_px, v[i, PXI])
		v[i, PYI] = vifelse(alive, new_py, v[i, PYI])

		e_factor = 0.25 / rel_p
		energy_change = e_factor * (kick_x * (kick_x + 2 * px_old) +
									kick_y * (kick_y + 2 * py_old)) + 
									0.5 * (dnk[1][1] * dsigma_ds[1] * sigma[1] + dnk[2][2] * dsigma_ds[2] * sigma[2])

		new_pz = v[i, PZI] + energy_change

		v[i, PZI] = vifelse(alive, new_pz, v[i, PZI])

		rel_p = 1 + v[i, PZI]
		pc = rel_p*p0c 
		#change to coord species later
		mc2 = BeamTracking.massof(Species("electron"))
		E_tot = sqrt(pc*pc + mc2*mc2)

		tilde_m = mc2 / p0c
		new_beta = pc / E_tot
		new_z = v[i, ZI] * (new_beta / beta_0)



		v[i, ZI] = vifelse(alive, new_z, v[i, ZI])
        

        #new stuff
        rel_p = 1 + v[i, PZI]
        ps_02 = rel_p * rel_p - v[i, PXI] * v[i, PXI] - v[i, PYI] * v[i, PYI]
        pc = rel_p*p0c 
        E_tot = sqrt(pc*pc + mc2*mc2)

        tilde_m = mc2 / p0c
        beta_0 = pc / E_tot
        gamsqr_0 = 1 / (1 - beta_0 * beta_0)
	end
	# v[i,:] = offset_particle(i, v, x_offset, y_offset, 
    #                     z_offset, x_pitch, y_pitch, 
    #                     tilt, false)
    print("v[i,XI] before last drift: ", v[i,XI], "\n")
	exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, -s_lab)
    print("v[i,XI] after beam-beam: ", v[i,XI], "\n")
end

"""
    strong_beam_center(z, crab, crab_tilt)

Calculate the (x,y,z) position of a slice within the strong beam in body coordinates.

# Arguments
- `z::Float64`: Position along the beam. Positive z is at the head of the bunch.
- `crab::Vector{Float64}`: Crab coefficients [crab_x1, crab_x2, crab_x3, crab_x4, crab_x5]
- `crab_tilt::Float64`: Crab tilt angle in radians

# Returns
- `center::Vector{Float64}`: 3-element vector containing (x, y, z) position of slice in body coordinates.
                             Positive z_strong is at the tail of the bunch.

# Notes
The function transforms z to z_strong = -z, so positive z_strong corresponds to the tail of the bunch.
"""
function strong_beam_center(z::Float64, crab::Vector{Float64}, crab_tilt::Float64)
    z_strong = -z
    
    # Using Horner's method for numerical stability
    r = (((crab[5] * z_strong + crab[4]) * z_strong + crab[3]) * z_strong + 
         crab[2]) * z_strong + crab[1]
    
    # Calculate center position
    center_x = z_strong * sin(r) * cos(crab_tilt)
    center_y = z_strong * sin(r) * sin(crab_tilt)
    center_z = -z_strong * cos(r)
    
    return [center_x, center_y, center_z]
end

# """
#     super_zbrent(func, x1, x2, rel_tol, abs_tol)

# Julia implementation of Brent's method for root finding, based on Bmad's super_zbrent.
# Finds the root of a function within the bracket [x1, x2].

# The x-tolerance is: x-tolerance = |x_root| * rel_tol + abs_tol

# ## Arguments
# - `func`:     Function whose root is to be found. func(x) should return (value, status)
# - `x1`, `x2`: Bracket values (x2 does not have to be larger than x1)
# - `rel_tol`:  Relative tolerance for root
# - `abs_tol`:  Absolute tolerance for root

# ## Returns
# - `(x_zero, status)` where status is:
#     - 0:  Normal convergence
#     - -1: Root not bracketed
#     - -2: Max iterations exceeded
#     - Other: Set by func
# """
# function super_zbrent(func, x1, x2, rel_tol, abs_tol)
#     itmax = 100
    
#     a = x1
#     print("a: ", a, "\n")
#     b = x2
#     print("b: ", b, "\n")
#     fa, status = func(a)
#     print("fa: ", fa, " status: ", status, "\n")
#     if status != 0
#         return NaN, status
#     end
#     fb, status = func(b)
#     print("fb: ", fb, " status: ", status, "\n")
#     if status != 0
#         return NaN, status
#     end
    
#     # Check if root is bracketed
#     if (fa > 0 && fb > 0) || (fa < 0 && fb < 0)
#         return NaN, -1
#     end
    
#     c = b
#     print("c: ", c, "\n")
#     fc = fb
#     print("fc: ", fc, "\n")
#     d = 0.0
#     e = 0.0
    
#     for iter in 1:itmax
#         if (fb > 0 && fc > 0) || (fb < 0 && fc < 0)
#             c = a
#             print("c updated: ", c, "\n")
#             fc = fa
#             print("fc updated: ", fc, "\n")
#             d = b - a
#             print("d updated: ", d, "\n")
#             e = d
#             print("e updated: ", e, "\n")
#         end
        
#         if abs(fc) < abs(fb)
#             a = b
#             print("a updated: ", a, "\n")
#             b = c
#             print("b updated: ", b, "\n")
#             c = a
#             print("c updated: ", c, "\n")
#             fa = fb
#             print("fa updated: ", fa, "\n")
#             fb = fc
#             print("fb updated: ", fb, "\n")
#             fc = fa
#             print("fc updated: ", fc, "\n")
#         end
        
#         tol1 = 0.5 * (rel_tol * abs(b) + abs_tol)
#         print("tol1: ", tol1, "\n")
#         xm = 0.5 * (c - b)
#         print("xm: ", xm, "\n")
        
#         if abs(xm) <= tol1 || fb == 0.0
#             # Converged
#             if fb == 0
#                 print("Root found exactly at b: ", b, "\n")
#                 return b, 0
#             else
#                 # Linear interpolation for final value
#                 print("Converged with tolerance at b: ", b, "\n")
#                 x_zero = (b * fc - c * fb) / (fc - fb)
#                 return x_zero, 0
#             end
#         end
        
#         if abs(e) >= tol1 && abs(fa) > abs(fb)
#             s = fb / fa
#             print("s: ", s, "\n")
#             if a == c
#                 p = 2.0 * xm * s
#                 print("p (secant): ", p, "\n")
#                 q = 1.0 - s
#                 print("q (secant): ", q, "\n")
#             else
#                 q = fa / fc
#                 print("q (inverse quadratic): ", q, "\n")
#                 r = fb / fc
#                 print("r (inverse quadratic): ", r, "\n")
#                 p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
#                 print("p (inverse quadratic): ", p, "\n")
#                 q = (q - 1.0) * (r - 1.0) * (s - 1.0)
#                 print("q (inverse quadratic updated): ", q, "\n")
#             end
            
#             if p > 0.0
#                 q = -q
#                 print("q adjusted: ", q, "\n")
#             end
#             p = abs(p)
#             print("p abs: ", p, "\n")
            
#             if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))
#                 e = d
#                 print("e updated for interpolation: ", e, "\n")
#                 d = p / q
#                 print("d updated for interpolation: ", d, "\n")
#             else
#                 d = xm
#                 print("d set to xm: ", d, "\n")
#                 e = d
#                 print("e set to d: ", e, "\n")
#             end
#         else
#             d = xm
#             print("d set to xm (bisection): ", d, "\n")
#             e = d
#             print("e set to d (bisection): ", e, "\n")
#         end
        
#         a = b
#         print("a updated to b: ", a, "\n")
#         fa = fb
#         print("fa updated to fb: ", fa, "\n")
        
#         if abs(d) > tol1
#             b = b + d
#             print("b updated with d: ", b, "\n")
#         else
#             b = b + sign(xm) * tol1
#             print("b updated with tol1: ", b, "\n")
#         end
        
#         fb, status = func(b)
#         print("new fb: ", fb, " status: ", status, "\n")
#         if status != 0
#             print("Function status error: ", status, "\n")
#             return b, status
#         end
#     end
    
#     # Max iterations exceeded
#     print("Max iterations exceeded\n")
#     return b, -2
# end

"""
    super_zbrent(func, x1, x2, rel_tol, abs_tol)

SIMD-vectorized Brent's method for root finding, based on Bmad's super_zbrent.
Finds roots for multiple particles simultaneously, continuing until all particles converge.

The x-tolerance is: x-tolerance = |x_root| * rel_tol + abs_tol

## Arguments
- `func`:     Function whose root is to be found. func(x) should return (value, status)
              Both x and value can be SIMD vectors
- `x1`, `x2`: Bracket values (can be scalars or SIMD vectors)
- `rel_tol`:  Relative tolerance for root
- `abs_tol`:  Absolute tolerance for root

## Returns
- `(x_zero, status)` where status is:
    - 0:  Normal convergence (all particles)
    - -1: Root not bracketed (any particle)
    - -2: Max iterations exceeded (any particle)
    - Other: Set by func
"""
function super_zbrent(func, x1, x2, rel_tol, abs_tol)
    itmax = 100
    
    a = x1
    b = x2
    fa, status = func(a)
    if status != 0
        return typeof(a)(NaN), status
    end
    fb, status = func(b)
    if status != 0
        return typeof(b)(NaN), status
    end
    z0 = zero(fa) # Vec of zeros of same type
    # Check if root is bracketed (for each particle in SIMD vector)
    not_bracketed = ((fa > z0) & (fb > z0)) | ((fa < z0) & (fb < z0))
    if any(not_bracketed)
        return typeof(b)(NaN), -1
    end
    
    c = b
    fc = fb
    d = zero(b)
    e = zero(b)
    
    # Track convergence for each particle (as SIMD Vec)
    converged = zero(b) == zero(b)  # Creates Vec{N, Bool} with all false
    
    for iter in 1:itmax
        # Update c when signs match
        same_sign = ((fb > z0) & (fc > z0)) | ((fb < z0) & (fc < z0))
        c = vifelse(same_sign, a, c)
        fc = vifelse(same_sign, fa, fc)
        d = vifelse(same_sign, b - a, d)
        e = vifelse(same_sign, d, e)
        
        # Swap if c is better estimate than b
        swap = abs(fc) < abs(fb)
        a_new = vifelse(swap, b, a)
        b_new = vifelse(swap, c, b)
        c_new = vifelse(swap, a_new, c)
        fa_new = vifelse(swap, fb, fa)
        fb_new = vifelse(swap, fc, fb)
        fc_new = vifelse(swap, fa_new, fc)
        
        a = a_new
        b = b_new
        c = c_new
        fa = fa_new
        fb = fb_new
        fc = fc_new
        
        # Convergence check
        tol1 = 0.5 * (rel_tol * abs(b) + abs_tol)
        print("tol1: ", tol1, "\n")
        xm = 0.5 * (c - b)
        print("xm: ", xm, "\n")
        
        # Check which particles have converged
        converged_this_iter = (abs(xm) <= tol1) | (fb == z0)
        converged = converged | converged_this_iter
        
        # If all particles converged, return
        if all(converged)
            # Linear interpolation for final value
            x_zero = (b * fc - c * fb) / (fc - fb)
            x_zero = vifelse(fb == z0, b, x_zero)
            return x_zero, 0
        end
        
        # Interpolation or bisection
        use_interp = (abs(e) >= tol1) & (abs(fa) > abs(fb))
        
        s = fb / fa
        print("s: ", s, "\n")
        
        # Secant vs inverse quadratic interpolation
        is_secant = a == c
        
        # Secant method
        p_secant = 2.0 * xm * s
        q_secant = 1.0 - s
        
        # Inverse quadratic interpolation
        q_iq = fa / fc
        r_iq = fb / fc
        p_iq = s * (2.0 * xm * q_iq * (q_iq - r_iq) - (b - a) * (r_iq - 1.0))
        q_iq_final = (q_iq - 1.0) * (r_iq - 1.0) * (s - 1.0)
        
        p = vifelse(is_secant, p_secant, p_iq)
        q = vifelse(is_secant, q_secant, q_iq_final)
        
        # Adjust q sign
        q = vifelse(p > 0.0, -q, q)
        p = abs(p)
        
        # Check if interpolation is acceptable
        accept_interp = 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))
        
        e_new = vifelse(use_interp & accept_interp, d, xm)
        d_new = vifelse(use_interp & accept_interp, p / q, xm)
        
        e = e_new
        d = d_new
        
        # Move to next point
        a = b
        fa = fb
        
        # Update b
        large_step = abs(d) > tol1
        b = vifelse(large_step, b + d, b + sign(xm) * tol1)
        
        fb, status = func(b)
        if status != 0
            print("Function status error: ", status, "\n")
            return b, status
        end
    end
    
    # Max iterations exceeded
    print("Max iterations exceeded\n")
    return b, -2
end

"""
    find_collision_point(i, coords, v, slice_center, beta_0, gamsqr_0, tilde_m, 
                        beta_strong, s_lab_current, s0_factor)

Find the s-position where the weak particle and strong beam slice meet using Brent's method.

## Arguments
- `i`:             Particle index
- `coords`:        Coordinate structure
- `v`:             Particle coordinate array
- `slice_center`:  [x,y,z] position of strong beam slice
- `beta_0`:        Reference beta
- `gamsqr_0`:      Squared Lorentz factor
- `tilde_m`:       Normalized mass
- `beta_strong`:   Strong beam beta
- `s_lab_current`: Current s-position in lab frame
- `s0_factor`:     Velocity ratio factor

## Returns
- `s_lab`:         Lab frame s-position of collision
"""
function find_collision_point(i, coords, v, slice_center, 
                             beta_0, gamsqr_0, tilde_m, beta_strong,
                             s_lab_current, s0_factor,s00, z_offset, part_time1, p0c,time0,part_time0)

	# Save initial state
	beta_l = ((p0c)*(sqrt((1 + v[i,PZI])*(1 + v[i,PZI]) - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI])))/sqrt((p0c*tilde_m)^2 + (p0c*(1 + v[i,PZI]))*(p0c*(1 + v[i,PZI])))
	p_l = sqrt((1 + v[i,PZI])*(1 + v[i,PZI]) - v[i,PXI]*v[i,PXI] - v[i,PYI]*v[i,PYI])
    p_rel = 1 + v[i,PZI]
    px_rel = v[i,PXI]/p_rel
    py_rel = v[i,PYI]/p_rel
    pxy2 = px_rel*px_rel + py_rel*py_rel
    # dt = length / (orb%beta * ps_rel * c_light)
	#beta_ref = 1/sqrt(1 + tilde_m*tilde_m);
	v_save = SA[v[i,XI], v[i,PXI], v[i,YI], v[i,PYI], v[i,ZI], v[i,PZI]]
	s_body0 = -z_offset + s_lab_current
	v[i,XI] = v_save[1]
	v[i,PXI] = v_save[2]
	v[i,YI] = v_save[3]
	v[i,PYI] = v_save[4]
	v[i,ZI] = v_save[5]
	v[i,PZI] = v_save[6]
	beta_0s = beta_0
	if typeof(v[i,XI]) <: TPS
		v[i,XI] = Float64(scalar.(v_save[1]))
		v[i,PXI] = Float64(scalar.(v_save[2]))
		v[i,YI] = Float64(scalar.(v_save[3]))
		v[i,PYI] = Float64(scalar.(v_save[4]))
		v[i,ZI] = Float64(scalar.(v_save[5]))
		v[i,PZI] = Float64(scalar.(v_save[6]))
		beta_0s = Float64(scalar(beta_0))
		p_l = Float64(scalar(p_l))
		s_body = Ref(Float64(scalar(s_body0)))
        s_lab = Ref(Float64(scalar(s_lab_current)))
		gamsqr_0s = Float64(scalar(gamsqr_0))
        ps_rel = Float64(scalar.(sqrt(1 - pxy2)))
        time = Ref(Float64(scalar(time0[])))
        part_time0 = Float64(scalar(part_time0))
	else
        time = Ref(time0[])
        ps_rel = sqrt(1 - pxy2)
		s_body = Ref(s_body0)
        s_lab = Ref(s_lab_current)
		gamsqr_0s = gamsqr_0
	end

	function collision_func(s_lab_target)
        print("\n")
		print("s_body before drift: ", s_body[], "\n")
        print("s_lab: ", s_lab[], "\n")
		del_s = s_lab_target - s_lab[]
        print("s_lab_target: ", s_lab_target, "\n")
        print("del_s: ", del_s, "\n")
		exact_drift!(i, coords, beta_0s, gamsqr_0, tilde_m, del_s)
        dt = del_s / (beta_0s * C_LIGHT* ps_rel)
        print("time before dt: ", time[] + part_time0, "\n")
        print("dt: ", dt, "\n")
        time[] = time[] + dt
        print("part_time0: ", part_time0, "\n")
        print("time: ", time[], "\n")
        s_lab[] = s_lab_target
		s_body[] = s_body[] + del_s
        print("s_body after drift: ", s_body[], "\n")
        
		s_slice = slice_center[3] - beta_strong * C_LIGHT * (time[])
		print("s_slice: ", s_slice, "\n")
		return s_body[] - s_slice, 0
	end
	# Initial guess for s_lab
	if typeof(v[i,XI]) <: TPS
		s0 = Float64(scalar(s00 + slice_center[3] * s0_factor))
	
    else
        s0 = s00 + slice_center[3] * s0_factor
    end
	ds = max(abs(slice_center[3]) * 0.1, 1.0) 
	
	s_lab, status = super_zbrent(collision_func, s0 - ds, s0 + ds, 1e-12, 1e-12)
    print("s_lab: ", s_lab, " status: ", status, "\n") 
	if status != 0
		# If root finding failed, return to initial guess
		s_lab = s0
	end
	if typeof(v[i,XI]) <: TPS
        dt = (s_lab[] - s_lab_current) / (beta_0 * C_LIGHT * sqrt(1 - pxy2))
		s_lab = slice_center[3] - beta_strong * C_LIGHT * (time0[] + dt)
	end

	# Restore original coordinates regardless of status
	v[i,XI] = v_save[1]
	v[i,PXI] = v_save[2]
	v[i,YI] = v_save[3]
	v[i,PYI] = v_save[4]
	v[i,ZI] = v_save[5]
	v[i,PZI] = v_save[6]

	return s_lab
	
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

	z_positions = zeros(n_slices)
	for i ∈ 1:n_slices
		prob = (i - 0.5) / n_slices - 0.5

        # erf returns a slightly different number than Fortran ERF, idk why
        z_norm = inverse( x -> erf(x/sqrt(2.0))/2.0, prob, -5.0, 5.0, 1.0e-5)
		z_positions[i] = sig_z * z_norm
	end
    print("z_positions: ", z_positions, "\n")
	return z_positions
end
"""
    inverse(funct, y, x1, x2, tol)

Find the inverse of a function Y = funct(X) to return X given Y using Brent's method.

This is a root-finding algorithm that combines bisection, secant, and inverse 
quadratic interpolation methods. It is a slight modification of the ZBRENT 
function from Numerical Recipes.

# Arguments
- `funct::Function`: Function to be inverted (must take a single Float64 argument)
- `y::Float64`: Value to be inverted (find x such that funct(x) = y)
- `x1::Float64`: Lower bound for search (x must satisfy x1 < x < x2)
- `x2::Float64`: Upper bound for search
- `tol::Float64`: Absolute accuracy tolerance for x

# Returns
- `x::Float64`: Solution such that funct(x) ≈ y

# Throws
- `ErrorException`: If root is not bracketed (funct(x1) and funct(x2) must have opposite signs relative to y)
- `@warn`: If maximum iterations exceeded (returns best estimate)

# Examples
```julia
# Find x such that x^2 = 4 (should return 2.0)
f(x) = x^2
x = inverse(f, 4.0, 0.0, 5.0, 1e-6)

# Invert exponential function
x = inverse(exp, 7.389, 0.0, 5.0, 1e-8)
```
"""
function inverse(funct::Function, y::Float64, x1::Float64, x2::Float64, tol::Float64)
    itmax = 100
    eps = 3.0e-8
    
    a = x1
    b = x2
    fa = funct(a) - y
    fb = funct(b) - y
    
    # Check that root is bracketed
    if fb * fa > 0.0
        error("Root must be bracketed for inverse function. " *
              "funct(x1)-y and funct(x2)-y must have opposite signs.")
    end
    
    c = b
    fc = fb
    d = 0.0
    e = 0.0
    
    for iter in 1:itmax
        # Ensure that b is the best estimate
        if fb * fc > 0.0
            c = a
            fc = fa
            d = b - a
            e = d
        end
        
        # If c is a better estimate than b, swap them
        if abs(fc) < abs(fb)
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        
        # Convergence check
        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        
        if abs(xm) <= tol1 || fb == 0.0
            return b
        end
        
        # Attempt inverse quadratic interpolation or secant method
        if abs(e) >= tol1 && abs(fa) > abs(fb)
            s = fb / fa
            
            if a == c
                # Secant method
                p = 2.0 * xm * s
                q = 1.0 - s
            else
                # Inverse quadratic interpolation
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            end
            
            # Ensure p is positive
            if p > 0.0
                q = -q
            end
            p = abs(p)
            
            # Check if interpolation is acceptable
            if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))
                # Accept interpolation
                e = d
                d = p / q
            else
                # Use bisection
                d = xm
                e = d
            end
        else
            # Use bisection
            d = xm
            e = d
        end
        
        # Move to next point
        a = b
        fa = fb
        
        if abs(d) > tol1
            b = b + d
        else
            b = b + sign(xm)*abs(tol1)
        end
        
        fb = funct(b) - y
    end
    
    @warn "inverse: Maximum iterations ($itmax) exceeded. Returning best estimate."
    return b
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
    #function to compute (1 - exp(-a)) / a using Taylor expansion for small a
	function one_minus_exp_over_a(a; nterms=10)
		s = zero(a)
		ak = one(a)          # a^0
		for k in 0:(nterms-1)
			s += ((-1)^k / factorial(k+1)) * ak
			ak *= a 
		end
		return s
	end
    sig_x, sig_y = sigma
    print("simgma: ", sigma, "\n")

    r_orig     = sig_y / sig_x
    x_norm_orig = dx / sig_x
    y_norm_orig = dy / sig_y

    # --- Round beam mask ---
    round_mask = abs(r_orig - 1.0) < 0.001
    amp = (x_norm_orig*x_norm_orig + y_norm_orig*y_norm_orig) / 2
    if typeof(amp) <: TPS
        scale_round = 4 * pi * one_minus_exp_over_a(amp, nterms = amp.mo)
    else
        scale_round = 4 * pi * (1 - exp(-amp)) / amp
        scale_round = vifelse(amp > 30, 4 * pi / amp, scale_round)
        scale_round = vifelse(amp < 1e-4, 4 * pi, scale_round)
    end

    # Round beam branch
    nkx_round = -x_norm_orig * scale_round
    nky_round = -y_norm_orig * scale_round

    # --- Flipped / elliptical branch ---
    flipped   = r_orig > 1.0
    r_f       = vifelse(flipped, 1.0 / r_orig, r_orig)
    x_norm_f  = vifelse(flipped, y_norm_orig, x_norm_orig)
    y_norm_f  = vifelse(flipped, x_norm_orig, y_norm_orig)
    sig_x_f   = vifelse(flipped, sig_y, sig_x)
    print("sig_x_f: ", sig_x_f, "\n")
    sig_y_f   = vifelse(flipped, sig_x, sig_y)
    print("sig_y_f: ", sig_y_f, "\n")

    denom = vifelse(round_mask, 1.0, 1.0 / sqrt(2 * (1 - r_f*r_f) + 1e-12))
    print("denom: ", denom, "\n")
    scale_prefactor = 4 * sqrt(pi^3) * denom * (1 + r_f)
    print("scale_prefactor: ", scale_prefactor, "\n")

    u  = abs(x_norm_f) * denom
    v  = abs(y_norm_f) * denom
    sx = sign(x_norm_f)
    print("sx: ", sx, "\n")
    sy = sign(y_norm_f)
    print("sy: ", sy, "\n")

    # --- Split-complex representation ---
    zr1 = u
    zi1 = r_f * v
    zr2 = r_f * u
    zi2 = v

    N = length(dx)

    # Use NTuple for lane-wise Faddeeva results
    f_real = ntuple(_ -> 0.0, N)
    f_imag = ntuple(_ -> 0.0, N)
    dw1x_dx = ntuple(_ -> 0.0, N)
    dw1y_dx = ntuple(_ -> 0.0, N)
    dw1x_dy = ntuple(_ -> 0.0, N)
    dw1y_dy = ntuple(_ -> 0.0, N)
    dw2x_dx = ntuple(_ -> 0.0, N)
    dw2y_dx = ntuple(_ -> 0.0, N)
    dw2x_dy = ntuple(_ -> 0.0, N)
    dw2y_dy = ntuple(_ -> 0.0, N)
    w1real = ntuple(_ -> 0.0, N)
    w1imag = ntuple(_ -> 0.0, N)
    w2real = ntuple(_ -> 0.0, N)
    w2imag = ntuple(_ -> 0.0, N)
    df_du_real = ntuple(_ -> 0.0, N)
    df_du_imag = ntuple(_ -> 0.0, N)
    df_dv_real = ntuple(_ -> 0.0, N)
    df_dv_imag = ntuple(_ -> 0.0, N)
    # Compute Faddeeva per lane
    for i in 1:N
        z1 = Complex(zr1[i], zi1[i])
        z2 = Complex(zr2[i], zi2[i])
        w1 = exp(-z1^2) * erfc(-im*z1)
        print("w1: ", w1, "\n")
        w2 = exp(-z2^2) * erfc(-im*z2)
        print("w2: ", w2, "\n")
        arg = (1 - r_f[i]^2) * (u[i]^2 + v[i]^2)
        expon = exp(-arg)
        print("expon: ", expon, "\n")
        f = w1 - expon * w2
        print("f: ", f, "\n")
        f_real = Base.setindex(f_real, real(f), i)
        f_imag = Base.setindex(f_imag, imag(f), i)
        dw1_dz1 = -2.0 * z1 * w1 + 2.0im / sqrt(pi)
        print("dw1_dz1: ", dw1_dz1, "\n")
        dw2_dz2 = -2.0 * z2 * w2 + 2.0im / sqrt(pi)
        print("dw2_dz2: ", dw2_dz2, "\n")

        df_du = dw1_dz1 - expon * (r_f[i] * dw2_dz2 - 2.0 * (1 - r_f[i]^2) * u[i] * w2)
        print("df_du: ", df_du, "\n")
        df_dv = im * r_f[i] * dw1_dz1 - expon * (im * dw2_dz2 - 2.0 * (1 - r_f[i]^2) * v[i] * w2)
        print("df_dv: ", df_dv, "\n")
        df_du_real = Base.setindex(df_du_real, real(df_du), i)
        df_du_imag = Base.setindex(df_du_imag, imag(df_du), i)
        df_dv_real = Base.setindex(df_dv_real, real(df_dv), i)
        df_dv_imag = Base.setindex(df_dv_imag, imag(df_dv), i)
        dw1x_dx = Base.setindex(dw1x_dx, real(dw1_dz1), i)
        dw1y_dx = Base.setindex(dw1y_dx, imag(dw1_dz1), i)
        dw1x_dy = Base.setindex(dw1x_dy, -imag(dw1_dz1), i)
        dw1y_dy = Base.setindex(dw1y_dy, real(dw1_dz1), i)
        dw2x_dx = Base.setindex(dw2x_dx, real(dw2_dz2), i)
        dw2y_dx = Base.setindex(dw2y_dx, imag(dw2_dz2), i)
        dw2x_dy = Base.setindex(dw2x_dy, -imag(dw2_dz2), i)
        dw2y_dy = Base.setindex(dw2y_dy, real(dw2_dz2), i)
        w1real = Base.setindex(w1real, real(w1), i)
        w1imag = Base.setindex(w1imag, imag(w1), i)
        w2real = Base.setindex(w2real, real(w2), i)
        w2imag = Base.setindex(w2imag, imag(w2), i)

    end

    # Convert back to Vec
    if length(dx) > 1
        f_real_vec = SIMD.Vec(f_real...)
        f_imag_vec = SIMD.Vec(f_imag...)
        df_du_imag_vec = SIMD.Vec(df_du_imag...)
        df_dv_real_vec = SIMD.Vec(df_dv_real...)
        df_du_real_vec = SIMD.Vec(df_du_real...)
        df_dv_imag_vec = SIMD.Vec(df_dv_imag...)
        dw1x_dx_vec = SIMD.Vec(dw1x_dx...)
        dw1y_dx_vec = SIMD.Vec(dw1y_dx...)
        dw1x_dy_vec = SIMD.Vec(dw1x_dy...)
        dw1y_dy_vec = SIMD.Vec(dw1y_dy...)
        dw2x_dx_vec = SIMD.Vec(dw2x_dx...)
        dw2y_dx_vec = SIMD.Vec(dw2y_dx...)
        dw2x_dy_vec = SIMD.Vec(dw2x_dy...)
        dw2y_dy_vec = SIMD.Vec(dw2y_dy...)
        w1real_vec = SIMD.Vec(w1real...)
        w1imag_vec = SIMD.Vec(w1imag...)
        w2real_vec = SIMD.Vec(w2real...)
        w2imag_vec = SIMD.Vec(w2imag...)
    else
        f_real_vec = f_real[1]
        f_imag_vec = f_imag[1]
        df_du_imag_vec = df_du_imag[1]
        df_dv_real_vec = df_dv_real[1]
        df_du_real_vec = df_du_real[1]
        df_dv_imag_vec = df_dv_imag[1]
        dw1x_dx_vec = dw1x_dx[1]
        print("dw1x_dx_vec: ", dw1x_dx_vec, "\n")
        dw1y_dx_vec = dw1y_dx[1]
        print("dw1y_dx_vec: ", dw1y_dx_vec, "\n")
        dw1x_dy_vec = dw1x_dy[1]
        print("dw1x_dy_vec: ", dw1x_dy_vec, "\n")
        dw1y_dy_vec = dw1y_dy[1]
        print("dw1y_dy_vec: ", dw1y_dy_vec, "\n")
        dw2x_dx_vec = dw2x_dx[1]
        print("dw2x_dx_vec: ", dw2x_dx_vec, "\n")
        dw2y_dx_vec = dw2y_dx[1]
        print("dw2y_dx_vec: ", dw2y_dx_vec, "\n")
        dw2x_dy_vec = dw2x_dy[1]
        print("dw2x_dy_vec: ", dw2x_dy_vec, "\n")
        dw2y_dy_vec = dw2y_dy[1]
        print("dw2y_dy_vec: ", dw2y_dy_vec, "\n")
        w1real_vec = w1real[1]
        print("w1real_vec: ", w1real_vec, "\n")
        w1imag_vec = w1imag[1]
        print("w1imag_vec: ", w1imag_vec, "\n")
        w2real_vec = w2real[1]
        print("w2real_vec: ", w2real_vec, "\n")
        w2imag_vec = w2imag[1]
        print("w2imag_vec: ", w2imag_vec, "\n")
    end

    nkx_general = -scale_prefactor * sx * f_imag_vec
    nky_general = -scale_prefactor * sy * f_real_vec

    # --- Derivatives ---
    du_dx = sx * denom
    dv_dy = sy * denom

    # derivatives
    # dnk_x_dx = -scale_prefactor * sx * df_du_imag_vec * du_dx / sig_x_f
    # dnk_x_dy = -scale_prefactor * sx * df_dv_imag_vec * dv_dy / sig_y_f
    # dnk_y_dx = -scale_prefactor * sy * df_du_real_vec * du_dx / sig_x_f
    # dnk_y_dy = -scale_prefactor * sy * df_dv_real_vec * dv_dy / sig_y_f
    arg = (1 - r_f*r_f) * (u*u + v*v)
    expon = exp(-arg)
    dnk_x_dx = -scale_prefactor * denom * (dw1y_dx_vec/sig_x_f - expon * dw2y_dx_vec * r_f / sig_x_f + 2.0 * (1 - r_f*r_f) * w2imag_vec * expon * u/sig_x_f)
    dnk_x_dy = -scale_prefactor * denom * sx*sy* (dw1y_dy_vec*r_f/sig_y_f - expon * dw2y_dy_vec/sig_y_f + 2.0 * (1 - r_f*r_f) * w2imag_vec * expon * v/sig_y_f)
    dnk_y_dx = -scale_prefactor * denom * sx*sy*(dw1x_dx_vec/sig_x_f - expon * dw2x_dx_vec*r_f/sig_x_f + 2.0 * (1 - r_f*r_f) * w2real_vec * expon * u/sig_x_f)
    dnk_y_dy = -scale_prefactor * denom * (dw1x_dy_vec*r_f/sig_y_f - expon * dw2x_dy_vec/sig_y_f + 2.0 * (1 - r_f*r_f) * w2real_vec * expon * v/sig_y_f)

    # --- Flip Jacobian if needed ---
    nkx_after_flip = vifelse(flipped, nky_general, nkx_general)
    nky_after_flip = vifelse(flipped, nkx_general, nky_general)

    e11 = vifelse(flipped, dnk_y_dy, dnk_x_dx)
    e12 = vifelse(flipped, dnk_y_dx, dnk_x_dy)
    e21 = vifelse(flipped, dnk_x_dy, dnk_y_dx)
    e22 = vifelse(flipped, dnk_x_dx, dnk_y_dy)

    # --- Select round beam ---
    nk_x = vifelse(round_mask, nkx_round, nkx_after_flip)
    nk_y = vifelse(round_mask, nky_round, nky_after_flip)

    e11_final = vifelse(round_mask, -scale_round, e11)
    e12_final = vifelse(round_mask, zero(e11_final), e12)
    e21_final = vifelse(round_mask, zero(e11_final), e21)
    e22_final = vifelse(round_mask, -scale_round, e22)

    dnk_final = (
        (e11_final,  e12_final),
        (e21_final,  e22_final)
        )

    return nk_x, nk_y, dnk_final
end


"""
    strong_beam_sigma_calc(s_pos, sig_x0, sig_y0, s_twiss_ref, 
                          beta_a_strong, alpha_a_strong, beta_a, alpha_a_prev, alpha_a,
                          beta_b_strong, alpha_b_strong, beta_b, alpha_b_prev, alpha_b,
                          beam_strength, classical_radius_factor, p0c)

Calculate the strong beam sigmas, strong beam centroid offsets (due to crabbing), 
and the BBI force constant for a beambeam element.

# Arguments
- `sig_x0`: Reference sigma in x
- `sig_y0`: Reference sigma in y
- `ds`: Distance slice has traveled
- `beta_a_strong`: Strong beam beta in plane a (use 0 to calculate from lattice)
- `alpha_a_strong`: Strong beam alpha in plane a
- `beta_b_strong`: Strong beam beta in plane b (use 0 to calculate from lattice)
- `alpha_b_strong`: Strong beam alpha in plane b
- `beam_strength`: Strong beam strength parameter
- `p0c`: Reference momentum times c

# Returns
- `sigma`: Strong beam (x, y) sigmas
- `bbi_const`: BBI kick scale factor
- `dsigma_ds`: sig_x and sig_y longitudinal derivatives
"""
function strong_beam_sigma_calc(sig_x0, sig_y0, 
                                ds,
                                beta_a_strong, alpha_a_strong,
                                beta_b_strong, alpha_b_strong,
                                beam_strength, 
                                p0c)
    
    
    beta_a0 = beta_a_strong
    alpha_a0 = alpha_a_strong
    if beta_a0 == 0.0 
        sigma_x = sig_x0
        dsigma_ds_x = 0.0
    else
        gamma0 = (1.0 + alpha_a0^2) / beta_a0
        beta = beta_a0 - 2.0 * alpha_a0 * ds + gamma0 * ds^2
        sigma_x = sig_x0 * sqrt(beta / beta_a0)
        dsigma_ds_x = -(alpha_a0 - ds * gamma0) * sigma_x / beta
    end
    
    beta_b0 = beta_b_strong
    alpha_b0 = alpha_b_strong
    
    if beta_b0 == 0.0  
        sigma_y = sig_y0
        dsigma_ds_y = 0.0
    else
        gamma0 = (1.0 + alpha_b0^2) / beta_b0
        beta = beta_b0 - 2.0 * alpha_b0 * ds + gamma0 * ds^2
        sigma_y = sig_y0 * sqrt(beta / beta_b0)
        dsigma_ds_y = -(alpha_b0 - ds * gamma0) * sigma_y / beta
    end
    
    # Calculate BBI constant
    bbi_const = -beam_strength / 
                (2.0 * π * p0c * (sigma_x + sigma_y))
    
    return SA[sigma_x, sigma_y], bbi_const, SA[dsigma_ds_x, dsigma_ds_y]
end

"""
    offset_particle(orbit, x_offset, y_offset, z_offset, x_pitch, y_pitch, tilt, set::Bool)

Transform particle coordinates between laboratory and element body coordinates
accounting for misalignments and tilt.

# Arguments
- `orbit::Vector{Float64}`: 6-element particle coordinate vector [x, px, y, py, z, pz]
- `x_offset::Float64`: Horizontal offset
- `y_offset::Float64`: Vertical offset  
- `z_offset::Float64`: Longitudinal offset
- `x_pitch::Float64`: Pitch angle around x-axis (radians)
- `y_pitch::Float64`: Pitch angle around y-axis (radians)
- `tilt::Float64`: Tilt/roll angle around z-axis (radians)
- `set::Bool`: true = lab to body coords, false = body to lab coords

# Returns
- `orbit_out::Vector{Float64}`: Transformed 6-element coordinate vector

# Notes
- When set=true: Transforms from lab coordinates to element body coordinates
- When set=false: Transforms from element body coordinates back to lab coordinates
- The function applies offsets, pitches, and tilt in the appropriate order
"""
function offset_particle(i, orbit, x_offset::Float64, y_offset::Float64, 
                        z_offset::Float64, x_pitch::Float64, y_pitch::Float64, 
                        tilt::Float64, set::Bool)
    
    orbit_out = copy(orbit)
    x, px, y, py, z, pz = orbit_out[i, :]
    rel_p = 1.0 + pz
    
    if set  # Lab to body coordinates
        
        position_r = [x - x_offset, y - y_offset, z]
        
        # Create rotation matrix for pitch angles (inverse)
        ws = floor_angles_to_w_mat_inv(x_pitch, y_pitch, 0.0)
        
        # Rotate position
        position_r = ws * position_r
        
        pz_mag_sq = rel_p^2 - px^2 - py^2
        if pz_mag_sq <= 0
            error("Lost particle: pz^2 <= 0")
        end
        pz_mag = sqrt(pz_mag_sq)
        
        p_vec0 = [px, py, pz_mag]
        p_vec = ws * p_vec0
        
        orbit_out[1] = position_r[1]
        orbit_out[2] = p_vec[1]
        orbit_out[3] = position_r[2]
        orbit_out[4] = p_vec[2]
        orbit_out[5] = position_r[3]
        
        if tilt != 0.0
            orbit_out = tilt_coords(tilt, orbit_out)
        end
        
    else  # Body to lab coordinates
        
        # Unapply tilt first
        if tilt != 0.0
            orbit_out = tilt_coords(-tilt, orbit_out)
        end
        
        # Extract coordinates
        x, px, y, py, z, pz = orbit_out
        position_r = [x, y, z]
        
        # Create rotation matrix for pitch angles (forward)
        ws = floor_angles_to_w_mat(x_pitch, y_pitch, 0.0)
        
        # Calculate pz component
        pz_mag_sq = rel_p^2 - px^2 - py^2
        if pz_mag_sq <= 0
            error("Lost particle: pz^2 <= 0")
        end
        pz_mag = sqrt(pz_mag_sq)
        
        # Rotate momentum
        p_vec0 = [px, py, pz_mag]
        p_vec = ws * p_vec0
        
        # Rotate position and apply offsets
        position_r = ws * position_r
        position_r = position_r + [x_offset, y_offset, 0]
        
        # Update orbit
        orbit_out[1] = position_r[1]
        orbit_out[2] = p_vec[1]
        orbit_out[3] = position_r[2]
        orbit_out[4] = p_vec[2]
        orbit_out[5] = position_r[3]
    end
    
    return orbit_out
end

"""
    floor_angles_to_w_mat(x_pitch, y_pitch, z_pitch)

Create rotation matrix from floor angles (forward rotation).
"""
function floor_angles_to_w_mat(x_pitch::Float64, y_pitch::Float64, z_pitch::Float64)
    cx, sx = cos(x_pitch), sin(x_pitch)
    cy, sy = cos(y_pitch), sin(y_pitch)
    cz, sz = cos(z_pitch), sin(z_pitch)
    
    # Rotation matrix: Rz * Ry * Rx
    ws = [cy*cz                    -cy*sz                   sy;
          sx*sy*cz + cx*sz        -sx*sy*sz + cx*cz        -sx*cy;
          -cx*sy*cz + sx*sz        cx*sy*sz + sx*cz         cx*cy]
    
    return ws
end

"""
    floor_angles_to_w_mat_inv(x_pitch, y_pitch, z_pitch)

Create inverse rotation matrix from floor angles (inverse = transpose for rotation matrices).
"""
function floor_angles_to_w_mat_inv(x_pitch::Float64, y_pitch::Float64, z_pitch::Float64)
    return transpose(floor_angles_to_w_mat(x_pitch, y_pitch, z_pitch))
end

"""
    tilt_coords(tilt_angle, orbit)

Apply tilt rotation to particle coordinates.
"""
function tilt_coords(tilt_angle::Float64, orbit)
    orbit_out = copy(orbit)
    c = cos(tilt_angle)
    s = sin(tilt_angle)
    
    # Rotate position (x, y)
    x_new = c * orbit[1] + s * orbit[3]
    y_new = -s * orbit[1] + c * orbit[3]
    
    # Rotate momentum (px, py)
    px_new = c * orbit[2] + s * orbit[4]
    py_new = -s * orbit[2] + c * orbit[4]
    
    orbit_out[1] = x_new
    orbit_out[2] = px_new
    orbit_out[3] = y_new
    orbit_out[4] = py_new
    
    return orbit_out
end


# Deprecated version without SIMD
# function bbi_kick_faddeeva(dx, dy, sigma)
	
# 	#function to compute (1 - exp(-a)) / a using Taylor expansion for small a
# 	function one_minus_exp_over_a(a; nterms=10)
# 		s = zero(a)
# 		ak = one(a)          # a^0
# 		for k in 0:(nterms-1)
# 			s += ((-1)^k / factorial(k+1)) * ak
# 			ak *= a 
# 		end
# 		return s
# 	end
# 	sig_x, sig_y = sigma
#     print("sig_x: ", sig_x, " sig_y: ", sig_y, "\n")
# 	r = sig_y / sig_x
#     print("r: ", r, "\n")
# 	x_norm = dx / sig_x
#     print("x_norm: ", x_norm, "\n")
# 	y_norm = dy / sig_y
#     print("y_norm: ", y_norm, "\n")

# 	# Round beam case
# 	if abs(r - 1.0) < 0.001
# 		amp = (x_norm*x_norm + y_norm*y_norm) / 2

# 		if typeof(amp) <: TPS
# 			scale = 4 * pi * one_minus_exp_over_a(amp, nterms = amp.mo)
# 			dnk = SA[-scale 0.0;
# 		          0.0 -scale]
# 			return -x_norm * scale, -y_norm * scale, dnk
# 		end
# 		scale = 4 * pi * (1 - exp(-amp)) / amp
# 		scale = vifelse(amp > 30, 4 * pi / amp, scale)
# 		scale = vifelse(amp < 1e-4, 4 * pi, scale)
# 		dnk = SA[-scale 0.0;
# 		          0.0 -scale]
# 		return -x_norm * scale, -y_norm * scale, dnk
# 	end

# 	# Ensure r < 1 (swap if needed)
# 	flipped = (r > 1)
# 	if flipped
# 		r = 1 / r
# 		x_norm, y_norm = y_norm, x_norm
# 		sig_x, sig_y = sig_y, sig_x
# 	end

# 	denom = 1 / sqrt(2 * (1 - r^2))
# 	scale = 4 * sqrt(pi^3) * denom * (1 + r)
#     print("scale: ", scale, "\n")
# 	u = abs(x_norm) * denom
# 	v = abs(y_norm) * denom
# 	sx = sign(x_norm)
# 	sy = sign(y_norm)

# 	z1 = complex(u, r * v)
# 	z2 = complex(r * u, v)

# 	# Faddeeva function: w(z) = exp(-z²) * erfc(-i*z)
# 	w1 = exp(-z1^2) * erfc(-im * z1)
# 	w2 = exp(-z2^2) * erfc(-im * z2)

# 	arg = (1 - r^2) * (u^2 + v^2)
# 	expon = exp(-arg)


# 	f = w1 - expon * w2

# 	nk_x = -scale * sx * imag(f)
# 	nk_y = -scale * sy * real(f)

#     # Compute derivatives using Faddeeva derivative: dw/dz = -2z*w(z) + 2i/sqrt(π)
	
# 	dw1_dz1 = -2.0 * z1 * w1 + 2.0im / sqrt(pi)
# 	dw2_dz2 = -2.0 * z2 * w2 + 2.0im / sqrt(pi)
	
# 	df_du = dw1_dz1 - expon * (r * dw2_dz2 - 2.0 * (1 - r*r) * u * w2)
# 	df_dv = im * r * dw1_dz1 - expon * (im * dw2_dz2 - 2.0 * (1 - r*r) * v * w2)
	
# 	# Convert to derivatives with respect to x_norm and y_norm
# 	# u = |x_norm| / denom, v = |y_norm| / denom
# 	# du/dx_norm = sign(x_norm) / denom
# 	# dv/dy_norm = sign(y_norm) / denom
	
# 	du_dx = sx * denom
# 	dv_dy = sy * denom
	
# 	# nk_x = -scale * sx * imag(f)
# 	# nk_y = -scale * sy * real(f)
	
# 	dnk_x_dx_norm = -scale * sx * imag(df_du) * du_dx
# 	dnk_x_dy_norm = -scale * sx * imag(df_dv) * dv_dy
# 	dnk_y_dx_norm = -scale * sy * real(df_du) * du_dx
# 	dnk_y_dy_norm = -scale * sy * real(df_dv) * dv_dy
	
# 	# Convert to derivatives with respect to dx and dy (not normalized)
# 	# x_norm = dx / sig_x  =>  dx_norm/dx = 1/sig_x
# 	# y_norm = dy / sig_y  =>  dy_norm/dy = 1/sig_y
	
# 	dnk_x_dx = dnk_x_dx_norm / sig_x
# 	dnk_x_dy = dnk_x_dy_norm / sig_y
# 	dnk_y_dx = dnk_y_dx_norm / sig_x
# 	dnk_y_dy = dnk_y_dy_norm / sig_y
	
# 	if flipped
# 		nk_x, nk_y = nk_y, nk_x
# 		# Swap the Jacobian rows and columns
# 		dnk = SA[dnk_y_dy dnk_y_dx;
# 		         dnk_x_dy dnk_x_dx]
# 	else
# 		dnk = SA[dnk_x_dx dnk_x_dy;
# 		         dnk_y_dx dnk_y_dy]
# 	end

#     print("dnk: ", dnk, "\n")

# 	return nk_x, nk_y, dnk
# end