@testset "ExactTracking" begin
    @testset "Utility functions" begin
        dx_rot = 0.1
        dy_rot = -0.1
        dz_rot = 0.2
  
        W = [cos(dy_rot) 0 sin(dy_rot); 0 1 0; -sin(dy_rot) 0 cos(dy_rot)] *
          [1 0 0; 0 cos(dx_rot) sin(dx_rot); 0 -sin(dx_rot) cos(dx_rot)] *
          [cos(dz_rot) -sin(dz_rot) 0; sin(dz_rot) cos(dz_rot) 0; 0 0 1]
  
        # Test w_matrix function
        @test all(ExactTracking.w_matrix(dx_rot, dy_rot, dz_rot) .== W)
  
        Winv = [cos(dz_rot) sin(dz_rot) 0; -sin(dz_rot) cos(dz_rot) 0; 0 0 1] *
             [1 0 0; 0 cos(dx_rot) -sin(dx_rot); 0 sin(dx_rot) cos(dx_rot)] *
             [cos(dy_rot) 0 -sin(dy_rot); 0 1 0; sin(dy_rot) 0 cos(dy_rot)]
  
        # Test w_inv_matrix function
        @test all(ExactTracking.w_inv_matrix(dx_rot, dy_rot, dz_rot) .== Winv)
    end
  
    @testset "Kernels" begin
        function coeffs_approx_equal(v1, v2, ϵ)
            n = GTPSA.numcoefs(v1[1])
            for i in 1:6
                for j in 0:n-1
                    c1, c2 = v1[i][j], v2[i][j]
                    @test abs(c1 - c2) <= max(ϵ, ϵ * (abs(c1) + abs(c2)))
                end
            end
            return true
        end

        function patch(::Type{T}) where {T}
            dt = T(1e-9)
            dx = T(2)
            dy = T(3)
            dz = T(4)
            winv = ExactTracking.w_inv_matrix(T(5),T(6),T(7))
            return dt, dx, dy, dz, winv
        end

        # GTPSA patch map comparison
        include("bmad_maps/patch.jl")
        dt, dx, dy, dz, winv = patch(TPS64{d_z})
        L = winv[3,1]*dx + winv[3,2]*dy + winv[3,3]*dz
        v = transpose(@vars(d_z))
        v = ExactTracking.patch!(1, v, zeros(eltype(v),1,9), L, 10e6, 510998.95,dt,dx,dy,dz,winv)

        coeffs_approx_equal(v_z, v, 5e-8)

        # GTPSA drift map comparison
        include("bmad_maps/drift.jl")
        L = TPS64{d_z}(1)
        v = transpose(@vars(d_z))
        p0c = 10e6
        mc2 = 510998.95
        tilde_m = mc2/p0c
        beta_0 = p0c/sqrt(p0c^2 + mc2^2)
        gamsqr_0 = 1/(1 - beta_0^2)

        v = ExactTracking.exact_drift!(1, v, zeros(eltype(v),1,2), L, tilde_m, gamsqr_0, beta_0)

        coeffs_approx_equal(v_z, v, 1e-8)
        
        # GTPSA solenoid map comparison
        include("bmad_maps/solenoid.jl")
        L  = TPS64{d_z}(1)
        ks = TPS64{d_z}(2)
        v = transpose(@vars(d_z))
        v = ExactTracking.exact_solenoid!(1, v, zeros(eltype(v),1,8), L, ks, 10e6, 510998.95)

        coeffs_approx_equal(v_z, v, 1e-8)
    end 
    
end