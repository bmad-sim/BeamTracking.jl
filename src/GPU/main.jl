import Pkg
Pkg.activate("$(ENV["HOME"])/.julia/dev/BeamTracking")
using BeamTracking
using KernelAbstractions

using CUDA
using CUDA.CUDAKernels
if CUDA.functional()
  const backend = CUDABackend()
  CUDA.allowscalar(false)
end
include("element.jl")
include("MatrixKick_GPUext.jl")

using JACC
using Test

using Plots

function track_drift_kernel_cuda!(x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
      MatrixKick_GPUext.track_drift!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    end
    nothing
end

function track_quad_kernel_cuda!(x, y, z, px, py, pz, beta_ref, tilde_m, gamsqr_ref, L, Bn1, brho_)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
      MatrixKick_GPUext.track_quad!(i, x, y, z, px, py, pz, beta_ref, tilde_m, gamsqr_ref, L, Bn1, brho_)
    end
    nothing
end

function track_kernel_cuda!(x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        lat)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
      track!(i, x, y, z, px, py, pz,
          brho, beta_gamma_ref,
          lat)
    end
    nothing
end

function track_kernel_cuda!(x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        lat)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(x)
      track!(i, x, y, z, px, py, pz,
          brho, beta_gamma_ref,
          lat)
    end
    nothing
end

@kernel function track_kernel!(x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        lat)
    i = @index(Global, Linear)
    track!(i, x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        lat)
end

function track!(i, x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        lat)
    tilde_m    = 1 / beta_gamma_ref
    gamsqr_ref = 1 + beta_gamma_ref^2
    beta_ref   = beta_gamma_ref / sqrt(gamsqr_ref)

    for ele in lat
      L = ele.L
      if ele isa MatrixKick.Drift
          MatrixKick_GPUext.track_drift!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
      elseif ele isa MatrixKick.Quadrupole
          MatrixKick_GPUext.track_quad!(i, x, y, z, px, py, pz, beta_ref, tilde_m, gamsqr_ref, L, ele.Bn1, brho)
      end
    end

    return nothing
end

function track!(i, x, y, z, px, py, pz,
        brho, beta_gamma_ref,
        ele::TrackingLattice.Element)
    L = ele.L
    tilde_m    = 1 / beta_gamma_ref
    gamsqr_ref = 1 + beta_gamma_ref^2
    beta_ref   = beta_gamma_ref / sqrt(gamsqr_ref)

    if ele isa MatrixKick.Drift
        MatrixKick_GPUext.track_drift!(i, x, y, z, px, py, pz, L, beta_ref, tilde_m, gamsqr_ref)
    elseif ele isa MatrixKick.Quadrupole
        MatrixKick_GPUext.track_quad!(i, x, y, z, px, py, pz, beta_ref, tilde_m, gamsqr_ref, L, ele.Bn1, brho)
    end


    return nothing
end

function main()
    if CUDA.functional()
      CUDA.versioninfo()
    end

    NCells = 160
    NParticles = 1024 * 10 * 6

    Ld = 5.75
    drft = MatrixKick.Drift(Ld)
    drft_h = MatrixKick.Drift(Ld/2)
    @info "Created drift of length $(drft.L) m"

    Lq = 0.5
    Bn1 = 4.0
    qf = MatrixKick.Quadrupole(Lq,  Bn1)
    qd = MatrixKick.Quadrupole(Lq, -Bn1)
    @info "Created quads with gradient ±$(qf.Bn1) T/m, length $(qf.L)"

    # Create lattice on device
    fodo = TrackingLattice.Element[
      drft_h
      qf
      drft
      qd
      drft_h
    ]
    @info "Created FODO cell, type of object is $(typeof(fodo))"

    Lattice = repeat(fodo, NCells);
    @info "Created FODO Lattice, type of object is $(typeof(fodo)), length $(length(Lattice))"

    Lattice_d = JACC.array(Lattice)
    @info "Lattice is on GPU, $(typeof(Lattice_d))"

    xi  = ([ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ] .* ones(1024*10)')[:]
    pxi = ([ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ] .* ones(1024*10)')[:]
    yi  = ([ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ] .* ones(1024*10)')[:]
    pyi = ([ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ] .* ones(1024*10)')[:]
    zi  = ([ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ] .* ones(1024*10)')[:]
    pzi = ([ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ] .* ones(1024*10)')[:]

    e_minus = Species("electron")
    mec2 = massof(e_minus) # 0.51099895069 MeV
    ek1 =   18.e9;  # eV
    bg1 = sqrt(ek1 / mec2 * (ek1 / mec2 + 2))
    bunch_d = (species = e_minus, beta_gamma_ref = bg1,
               x = JACC.array(xi),
               px = JACC.array(pxi),
               y = JACC.array(yi),
               py = JACC.array(pyi),
               z = JACC.array(zi),
               pz = JACC.array(pzi))
    brho_ = brho(massof(bunch_d.species), bunch_d.beta_gamma_ref, chargeof(bunch_d.species))

    xt = deepcopy(xi)
    yt = deepcopy(yi)
    zt = deepcopy(zi)
    pxt = deepcopy(pxi)
    pyt = deepcopy(pyi)
    pzt = deepcopy(pzi)

    display("threads")
    gr()
    plx = scatter(xi, pxi, color = :red, label = "0 turns")
    ply = scatter(yi, pyi, color = :red, label = "0 turns")
    plz = scatter(zi, pzi, color = :red, label = "0 turns")
    @time begin
      Threads.@threads for i in 1:NParticles
          track!(i, xt, yt, zt, pxt, pyt, pzt,
          brho_, bunch_d.beta_gamma_ref,
          Lattice)
      end
    end
    scatter!(plx, xt, pxt, color = :blue, label = "1 turn")
    scatter!(ply, yt, pyt, color = :blue, label = "1 turn")
    scatter!(plz, yt, pyt, color = :blue, label = "1 turn")
    @time begin
      Threads.@threads for i in 1:NParticles
          track!(i, xt, yt, zt, pxt, pyt, pzt,
          brho_, bunch_d.beta_gamma_ref,
          Lattice)
      end
    end
    scatter!(plx, xt, pxt, color = :green, label = "2 turn")
    scatter!(ply, yt, pyt, color = :green, label = "2 turn")
    scatter!(plz, yt, pyt, color = :green, label = "2 turn")
    @time begin
      Threads.@threads for i in 1:NParticles
          track!(i, xt, yt, zt, pxt, pyt, pzt,
          brho_, bunch_d.beta_gamma_ref,
          Lattice)
      end
    end
    scatter!(plx, xt, pxt, color = :black, label = "3 turn")
    scatter!(ply, yt, pyt, color = :black, label = "3 turn")
    scatter!(plz, yt, pyt, color = :black, label = "3 turn")

    plot(plx, ply, plz, layout = (1,3), size = (1000, 400))
    savefig("phspace.pdf")

    display("JACC")
    @time JACC.parallel_for(NParticles, track!,
                      bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz,
                      brho_, bunch_d.beta_gamma_ref,
                      Lattice_d)
    @time JACC.parallel_for(NParticles, track!,
                      bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz,
                      brho_, bunch_d.beta_gamma_ref,
                      Lattice_d)
    @time JACC.parallel_for(NParticles, track!,
                      bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz,
                      brho_, bunch_d.beta_gamma_ref,
                      Lattice_d)

    @test Array(bunch_d.x)  ≈  xt (rtol=5.e-13)
    @test Array(bunch_d.y)  ≈  yt (rtol=5.e-13)
    @test Array(bunch_d.z)  ≈  zt (rtol=5.e-13)
    @test Array(bunch_d.px) ≈ pxt (rtol=5.e-13)
    @test Array(bunch_d.py) ≈ pyt (rtol=5.e-13)
    @test Array(bunch_d.pz) ≈ pzt (rtol=5.e-13)

    backend = CUDABackend()
    x_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    px_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    y_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    py_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    z_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    pz_ka = KernelAbstractions.allocate(backend, eltype(xi), size(xi))
    KernelAbstractions.copyto!(backend, x_ka, xi)
    KernelAbstractions.copyto!(backend, px_ka, pxi)
    KernelAbstractions.copyto!(backend, y_ka, yi)
    KernelAbstractions.copyto!(backend, py_ka, pyi)
    KernelAbstractions.copyto!(backend, z_ka, zi)
    KernelAbstractions.copyto!(backend, pz_ka, pzi)

    display("Kernel abstractions")
    kernel! = track_kernel!(backend, (256,))
    linear_input = [i for i in 1:NParticles]
    linear_input_d = KernelAbstractions.allocate(backend, eltype(linear_input), size(linear_input))
    KernelAbstractions.copyto!(backend, linear_input_d, linear_input)

    @time begin
      kernel!(x_ka, y_ka, z_ka, px_ka, py_ka, pz_ka, brho_, bg1, Lattice_d, ndrange = size(x_ka))
      KernelAbstractions.synchronize(backend)
    end
    @time begin
      kernel!(x_ka, y_ka, z_ka, px_ka, py_ka, pz_ka, brho_, bg1, Lattice_d, ndrange = size(x_ka))
      KernelAbstractions.synchronize(backend)
    end
    @time begin
      kernel!(x_ka, y_ka, z_ka, px_ka, py_ka, pz_ka, brho_, bg1, Lattice_d, ndrange = size(x_ka))
      KernelAbstractions.synchronize(backend)
    end

    @test Array(x_ka)  ≈  xt (rtol=5.e-13)
    @test Array(y_ka)  ≈  yt (rtol=5.e-13)
    @test Array(z_ka)  ≈  zt (rtol=5.e-13)
    @test Array(px_ka) ≈ pxt (rtol=5.e-13)
    @test Array(py_ka) ≈ pyt (rtol=5.e-13)
    @test Array(pz_ka) ≈ pzt (rtol=5.e-13)

    display("CUDA")

    bunch_d = (species = e_minus, beta_gamma_ref = bg1,
               x = CuArray(xi),
               px = CuArray(pxi),
               y = CuArray(yi),
               py = CuArray(pyi),
               z = CuArray(zi),
               pz = CuArray(pzi))
    @time begin
        CUDA.@sync @cuda threads=256 blocks=cld(length(xi), 256) track_kernel_cuda!(
           bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz, brho_, bg1, Lattice_d
          )
    end
    @time begin
        CUDA.@sync @cuda threads=256 blocks=cld(length(xi), 256) track_kernel_cuda!(
           bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz, brho_, bg1, Lattice_d
          )
    end
    @time begin
        CUDA.@sync @cuda threads=256 blocks=cld(length(xi), 256) track_kernel_cuda!(
           bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz, brho_, bg1, Lattice_d
          )
    end
    @test Array(bunch_d.x)  ≈  xt (rtol=5.e-13)
    @test Array(bunch_d.y)  ≈  yt (rtol=5.e-13)
    @test Array(bunch_d.z)  ≈  zt (rtol=5.e-13)
    @test Array(bunch_d.px) ≈ pxt (rtol=5.e-13)
    @test Array(bunch_d.py) ≈ pyt (rtol=5.e-13)
    @test Array(bunch_d.pz) ≈ pzt (rtol=5.e-13)

    display("MKLG(c)😂")
    bunch_d = (species = e_minus, beta_gamma_ref = bg1,
               x = CuArray(Float64, xi),
               px = CuArray(Float64, pxi),
               y = CuArray(Float64, yi),
               py = CuArray(Float64, pyi),
               z = CuArray(Float64, zi),
               pz = CuArray(Float64, pzi))
    @time begin
        tilde_m    = 1 / bg1
        gamsqr_ref = 1 + bg1^2
        beta_ref   = bg1 / sqrt(gamsqr_ref)
        graph = CUDA.capture(throw_error=true) do
          for ele in Lattice
            L = ele.L
            if ele isa MatrixKick.Drift
               @cuda threads=256 blocks=cld(length(xi), 256) track_drift_kernel_cuda!(bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz,
                                         L, beta_ref, tilde_m, gamsqr_ref)
            elseif ele isa MatrixKick.Quadrupole
               @cuda threads=256 blocks=cld(length(xi), 256) track_quad_kernel_cuda!(bunch_d.x, bunch_d.y, bunch_d.z, bunch_d.px, bunch_d.py, bunch_d.pz,
                                        beta_ref, tilde_m, gamsqr_ref, L, ele.Bn1, brho_)
            end
          end
       end
   end
   @time begin
    gr_ex = instantiate(graph)
   end
   @time begin
        CUDA.launch(gr_ex)
        CUDA.synchronize()
   end
    @time begin
        CUDA.launch(gr_ex)
        CUDA.synchronize()
   end
    @time begin
        CUDA.launch(gr_ex)
        CUDA.synchronize()
   end

   @test Array(bunch_d.x)  ≈  xt (rtol=5.e-13)
   @test Array(bunch_d.y)  ≈  yt (rtol=5.e-13)
   @test Array(bunch_d.z)  ≈  zt (rtol=5.e-13)
   @test Array(bunch_d.px) ≈ pxt (rtol=5.e-13)
   @test Array(bunch_d.py) ≈ pyt (rtol=5.e-13)
   @test Array(bunch_d.pz) ≈ pzt (rtol=5.e-13)
end

main()

