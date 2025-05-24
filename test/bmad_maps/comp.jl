using GTPSA
function nonzero_indices(v)
    [i-1 for i in 1:GTPSA.numcoefs(v) if abs(v[i-1]) > 2e-15]
end

function compare_indices(a, b)
set_a, set_b = Set(a), Set(b)
only_in_a = setdiff(set_a, set_b)
only_in_b = setdiff(set_b, set_a)
return only_in_a, only_in_b
end

function coeffs_approx_equal(v1, v2, ϵ)
n = min(GTPSA.numcoefs(v1), GTPSA.numcoefs(v2))
for i in 0:n-1
    c1, c2 = v1[i], v2[i]
    if c1 != 0 && c2 != 0
        if abs(c1 - c2) < ϵ * (abs(c1) + abs(c2))
            return true
        end
    end
end
return true
end

function index_to_monomial(t::TPS, idx::Integer, switch::Bool=false)
    # Handle scalar part (index 0)
    if idx == 0
        return zeros(Cuchar, GTPSA.numvars(t) + GTPSA.numparams(t))
    end
    
    # Create buffer for monomial orders
    mono = Vector{Cuchar}(undef, GTPSA.numvars(t) + GTPSA.numparams(t))
    
    # Create a reference for coefficient that cycle! will fill
    # (we don't actually use this value)
    coef = Ref{numtype(t)}()
    
    # For first-order monomials (common case), we can handle directly
    if idx <= GTPSA.numvars(t) + GTPSA.numparams(t)
        mono .= 0
        mono[idx] = 1
        if switch
            return true
        else
            return mono
        end
    end
    
    # For higher orders, use cycle! to find the monomial
    current_idx = cycle!(t, -1, length(mono), mono, coef)
    
    # Iterate until we reach our target index
    while current_idx >= 0 && current_idx < idx
        current_idx = cycle!(t, current_idx, length(mono), mono, coef)
    end
    
    # Check if we found the right index
    if current_idx == idx
        if switch
            return true
        else
            return coef.x, [Int64(x) for x in mono]
        end
    else
        if switch
            return false
        else
            error("Invalid monomial index: $idx, $current_idx")
        end
    end
end

using BeamTracking

include("patch.jl")
include("patch9.jl")
include("patch8.jl")

function patch(::Type{T}) where {T}
    dx = T(1)
    dy = T(2)
    dz = T(3)
    winv = ExactTracking.w_inv_matrix(T(4),T(5),T(6))
    return dx, dy, dz, winv
end

dx, dy, dz, winv = patch(Float64)

v = transpose(@vars(d_z9))
v = ExactTracking.patch!(1, v, zeros(eltype(v),1,9), 10e6,510998.95,dx,dy,dz,winv)

terms=Dict(1=>"x ",2=>"px",3=>"y ",4=>"py",5=>"z ",6=>"pz")

for i=1:6
    k_b = nonzero_indices(v_z9[i])
    k_s = nonzero_indices(v[i])
    B_terms,S_terms = compare_indices(k_b,k_s)
    println(terms[i]," coeffs only in Bmad/PTC:")
    for a in B_terms
        println(" ",index_to_monomial(v_z[i],a))
    end
    println()
    println()
    println(terms[i]," coeffs only in  SciBmad: ")
    for a in S_terms
        println(" ",index_to_monomial(v[i],a))
    end
    println()
    println()
    println("======================================================")
end