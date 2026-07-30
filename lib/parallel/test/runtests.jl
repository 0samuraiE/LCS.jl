using Test

using Adapt
using AMDGPU
using CUDA
using Parallel

backends = Tuple{String,KA.Backend}[("cpu", KA.CPU())]
AMDGPU.functional() && push!(backends, ("amdgpu", AMDGPU.ROCBackend()))
CUDA.has_cuda() && push!(backends, ("cuda", CUDA.CUDABackend()))

@testset "Parallel" begin
    @testset "foraxes" begin
        for (name, backend) in backends
            @testset "$name" begin
                let
                    A = KA.zeros(backend, Float64, 10)
                    Parallel.foraxes(backend, (1:10,)) do i
                        A[i] = Float64(i)
                    end
                    @test Array(A) == collect(1.0:10.0)
                end

                let
                    A = KA.zeros(backend, Int, 6, 8)
                    Parallel.foraxes(backend, (2:5, 3:7)) do i, j
                        A[i, j] = i + j
                    end
                    expected = zeros(Int, 6, 8)
                    for i in 2:5, j in 3:7
                        expected[i, j] = i + j
                    end
                    @test Array(A) == expected
                end

                let
                    A = KA.zeros(backend, Int, 4, 4, 4)
                    Parallel.foraxes(backend, (1:4, 1:4, 1:4)) do i, j, k
                        A[i, j, k] = i + j + k
                    end
                    expected = [i + j + k for i in 1:4, j in 1:4, k in 1:4]
                    @test Array(A) == expected
                end
            end
        end
    end

    @testset "coloring" begin
        for (name, backend) in backends
            @testset "$name" begin
                for axes in [(2:4, 2:4, 2:4), (1:4, 1:4, 1:4)]
                    @testset "axes=$(axes)" begin
                        ref = KA.zeros(backend, Int, 5, 5, 5)
                        Parallel.foraxes(backend, axes; coloring=Parallel.NoColoring()) do i, j, k
                            ref[i, j, k] = i * 100 + j * 10 + k
                        end

                        for coloring in [Parallel.RedBlack(), Parallel.RedBlackFast(), Parallel.RedBlackBlock(2)]
                            @testset "$(typeof(coloring))" begin
                                A = KA.zeros(backend, Int, 5, 5, 5)
                                Parallel.foraxes(backend, axes; coloring) do i, j, k
                                    A[i, j, k] = i * 100 + j * 10 + k
                                end
                                @test Array(A) == Array(ref)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "coloring reference" begin
        axes = (2:4, 2:4, 2:4)
        ref0x00 = [
            0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0
        ]
        ref0x01 = [
            0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 0 1 0 0; 0 1 0 1 0; 0 0 1 0 0; 0 0 0 0 0;;;
            0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0
        ]

        for (name, backend) in backends
            @testset "$name" begin
                for coloring in [Parallel.RedBlack(), Parallel.RedBlackFast()]
                    @testset "$(typeof(coloring))" begin
                        let
                            A = KA.zeros(backend, Int, 5, 5, 5)
                            Parallel._foraxes_colored(backend, axes, coloring, Val(0x00)) do i, j, k
                                A[i, j, k] = 1
                            end
                            @test Array(A) == ref0x00
                        end

                        let
                            A = KA.zeros(backend, Int, 5, 5, 5)
                            Parallel._foraxes_colored(backend, axes, coloring, Val(0x01)) do i, j, k
                                A[i, j, k] = 1
                            end
                            @test Array(A) == ref0x01
                        end
                    end
                end
            end
        end
    end

    @testset "hostaccess" begin
        for (name, backend) in backends
            @testset "$name" begin
                A = adapt(backend, collect(1:3))
                @test Parallel.hostaccess(A, 1) == 1
                @test Parallel.hostaccess(A, 1:1) == [1]
                @test Parallel.hostaccess(A, 1:2) == [1, 2]
                @test Parallel.hostaccess(A, 1:2:3) == [1, 3]
                @test Parallel.hostaccess(A, [1, 2, 3]) == [1, 2, 3]

                let
                    A = adapt(backend, reshape(collect(1:6), 2, 3))
                    @test Parallel.hostaccess(A, 1:2, 1:3) == [1 3 5; 2 4 6]
                    @test Parallel.hostaccess(A, 1:1, 2:3) == [3 5]
                end
            end
        end
    end
end
