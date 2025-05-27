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
        # GTPSA patch map comparison
        function patch_args(::Type{T}) where {T}
            p0c = T(10e6)
            mc2 = T(ELECTRON.mass)
            tilde_m = mc2/p0c
            gamsqr_0 = 1 + 1/tilde_m^2
            beta_0 = 1/sqrt(1 + tilde_m^2)
            dt = T(1e-9)
            dx = T(2)
            dy = T(3)
            dz = T(4)
            winv = ExactTracking.w_inv_matrix(T(5),T(6),T(7))
            L = winv[3,1]*dx + winv[3,2]*dy + winv[3,3]*dz
            return L, tilde_m, gamsqr_0, beta_0, dt, dx, dy, dz, winv
        end

        test_map(ExactTracking.patch!, "bmad_maps/patch.jl", patch_args(Float64)...; tol=5e-10)

        # GTPSA drift map comparison
        function drift_args(::Type{T}) where {T}
            L = T(1)
            p0c = T(10e6)
            mc2 = T(ELECTRON.mass)
            tilde_m = mc2/p0c
            gamsqr_0 = 1 + 1/tilde_m^2
            beta_0 = 1/sqrt(1 + tilde_m^2)
            return L, tilde_m, gamsqr_0, beta_0
        end

        test_map(ExactTracking.exact_drift!, "bmad_maps/drift.jl", drift_args(Float64)...; tol=5e-10)
        
        # GTPSA solenoid map comparison
        function solenoid_args(::Type{T}) where {T}
            L = T(1)
            ks = T(2)
            p0c = T(10e6)
            mc2 = T(ELECTRON.mass)
            tilde_m = mc2/p0c
            gamsqr_0 = 1 + 1/tilde_m^2
            beta_0 = 1/sqrt(1 + tilde_m^2)
            return L, ks, tilde_m, gamsqr_0, beta_0
        end

        test_map(ExactTracking.exact_solenoid!, "bmad_maps/solenoid.jl", solenoid_args(Float64)...; tol=5e-10)

    end 
    
end