using CUDA
Y = CuArray(zeros(2, 16 * 256))

function apa(Y)
    thread = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if thread <= size(Y)[2]
    @inbounds begin
        Y[1, thread] = randn()
        Y[2, thread] = rand() * 2π
    end
    end
    return
end

@cuda threads=16 blocks=256 apa(Y)

x = Array(Y[1,:]) .* cos.(Array(Y[2,:]))
y = Array(Y[1,:]) .* sin.(Array(Y[2,:]))

using Plots
scatter(x,y, markersize = 0.1)
savefig("res.png")