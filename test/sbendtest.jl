using Test
using BeamTracking.ExactTracking
using GTPSA, BeamTracking
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI

@kwdef struct SBend{T}
    dg :: T
    g :: T
    L :: T
    e1 :: T 
    e2 :: T 
end

@inline function exact_bend!(i, v, work, theta, gtot, g, L, mc2, p0c, beta_ref)
    @inbounds begin #@FastGTPSA! begin
        work[i,1] = sqrt((1 + v[i,PZI])^2 - v[i,PYI]^2) #pt
        work[i,2] = theta + asin(v[i,PXI] / work[i,1]) #phi1
        work[i,3] = gtot / work[i,1] #gp
        work[i,4] = 1+g*v[i,XI] # 1 + g * x
        work[i,5] = cos(work[i,2]) #cos(theta + phi1)
        work[i,6] = sincu(theta) #sincu(theta)
        work[i,7] = 2*work[i,4]*sin(work[i,2])*L*work[i,6]- work[i,3]*work[i,4]^2*L^2*(work[i,6])^2 #alpha
        if abs(work[i,2]) < π/2
            work[i,8] = work[i,7]/(sqrt(work[i,5]^2 + work[i,3]*work[i,7]) + work[i,5]) #xi
        else
            work[i,8] = (sqrt(work[i,5]^2 + work[i,3]*work[i,7]) - work[i,5]) / work[i,3] #xi
        end
        work[i,9] = -L*work[i,6]-v[i,XI]*sin(theta) #Lcv
        work[i,10] = 2 * (work[i,2] - atan(work[i,8], -work[i,9])) #theta_p
        work[i,11] = sqrt(work[i,9]^2 + work[i,8]^2) / sincu(work[i,10]/2) #Lp
        work[i,12] = (v[i,PZI] * p0c + p0c) #p
  
        v[i,XI] = v[i,XI]*cos(theta) - g/2*L^2*(sincu(theta/2))^2 + work[i,8]
        v[i,PXI] = work[i,1]*sin(work[i,2] - work[i,10])
        v[i,YI] = v[i,YI] + v[i,PYI]*work[i,11]/work[i,1] 
        v[i,ZI] = v[i,ZI] - (1 + v[i,PZI])*work[i,11]/work[i,1] + 
                        L*(work[i,12]/sqrt(mc2^2+work[i,12]^2))/beta_ref
    end #end
    return v
  end

function track!(bunch::Bunch, ele::SBend)
    i = 1
    work = zeros(typeof(BeamTracking.soaview(bunch)[i,1]), i, 12)
    theta = ele.g * ele.L
    gtot = ele.g + ele.dg 
    mc2 = BeamTracking.massof(bunch.species) * 10^6
    p0c = 10E6
    beta_ref = p0c / sqrt(mc2^2 + p0c^2)
    return exact_bend!(i, BeamTracking.soaview(bunch), work, theta, gtot, ele.g, ele.L, mc2, p0c, beta_ref)
end

xi  = 0.1
pxi = -7.5e-4
yi  = 0.1
pyi = -3e-4 	
zi  = 0.1
pzi = -1e-3
ps = [xi;; pxi;; yi;; pyi;; zi;; pzi]

d = Descriptor(6, 4)
Δx = @vars(d) 
ps = [Δx[1];; Δx[2];; Δx[3];; Δx[4];; Δx[5];; Δx[6]] + ps

b = BeamTracking.Bunch(ps, species = Species("electron"), Brho_ref=35000.)
ele = SBend(L = 0.5, g = 1.0, e1 = 0.0, e2 = 0.0, dg = 0.001)
x = track!(b, ele)
M = zeros(6,6)
for i = 1:6
    for j = 1:6
        M[i, j] = x[i][j]
    end
end

bmadstd = [  0.85374615     0.52189495     0.00000000     0.00004425     0.00000000     0.14734148
            -0.47990496     0.87794249     0.00000000     0.00014397     0.00000000     0.47942570   
            -0.00014415    -0.00003629     1.00000000     0.54837870     0.00000000     0.00015870
             0.00000000     0.00000000     0.00000000     1.00000000     0.00000000     0.00000000   
            -0.48001775    -0.12085250     0.00000000     0.00015870     1.00000000    -0.01861516 
             0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000]

@test M ≈ bmadstd