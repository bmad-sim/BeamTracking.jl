using LinearAlgebra

# Define a unit quaternion (45° rotation around z-axis)
θ = π / 4
q_r = cos(θ / 2)
q_i = 0.0
q_j = 0.0
q_k = sin(θ / 2)
q = [q_r, q_i, q_j, q_k]

function quat_to_rotmat(q)
    q_r, q_i, q_j, q_k = q
    R = [
        1-2*(q_j^2+q_k^2) 2*(q_i*q_j-q_k*q_r) 2*(q_i*q_k+q_j*q_r);
        2*(q_i*q_j+q_k*q_r) 1-2*(q_i^2+q_k^2) 2*(q_j*q_k-q_i*q_r);
        2*(q_i*q_k-q_j*q_r) 2*(q_j*q_k+q_i*q_r) 1-2*(q_i^2+q_j^2)
    ]
    return R
end

R = quat_to_rotmat(q)
eigvals(R)