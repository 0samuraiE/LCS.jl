using Test

using MPI
using Topologies

MPI.Initialized() || MPI.Init()

@testset "Topologies" begin
    @testset "single constructor" begin
        @test Topologies.Topology([0, 0, 0]; multiprocessing=false) ==
            Topologies.Topology([0, 0, 0]; multiprocessing=true)
    end

    topo = Topologies.Mock.mocktopology(; linear_rank=1, proc_dims=(3, 2, 2))

    @testset "proc_size" begin
        @test Topologies.proc_size(topo) == 3 * 2 * 2
    end

    @testset "proc_dims" begin
        @test Topologies.proc_dims(topo) == (3, 2, 2)
    end

    @testset "linear_rank" begin
        @test Topologies.linear_rank(topo) == 1
    end

    @testset "linear_ranks" begin
        expected = collect(0:(3 * 2 * 2 - 1))
        @test collect(Topologies.linear_ranks(topo)) == expected
    end

    @testset "cart_rank" begin
        @test Topologies.cart_rank(topo) == (0, 0, 1)
    end

    @testset "cart_to_linear_rank" begin
        @test Topologies.cart_to_linear_rank(topo, (0, 0, 0)) == 0
        @test Topologies.cart_to_linear_rank(topo, (2, 1, 1)) == 11
    end

    @testset "linear_to_cart_rank" begin
        @test Topologies.linear_to_cart_rank(topo, 0) == (0, 0, 0)
        @test Topologies.linear_to_cart_rank(topo, 1) == (0, 0, 1)
        @test Topologies.linear_to_cart_rank(topo, 11) == (2, 1, 1)
    end

    @testset "each_cart_rank" begin
        coords = collect(Topologies.each_cart_rank(topo))
        @test length(coords) == 12
        @test (0, 0, 0) in coords
        @test (2, 1, 1) in coords
    end

    @testset "isroot" begin
        @test !Topologies.isroot(topo)
        root = Topologies.Mock.mocktopology(; linear_rank=0, proc_dims=(3, 2, 2))
        @test Topologies.isroot(root)
    end

    @testset "processing" begin
        single = Topologies.Mock.mocktopology(; linear_rank=0, proc_dims=(1, 1, 1))
        @test Topologies.processing(single) isa Topologies.SingleProcessing
        @test !Topologies.is_multi_processing(single)
        @test Topologies.processing(topo) isa Topologies.MultiProcessing
        @test Topologies.is_multi_processing(topo)
    end

    @testset "mpi" begin
        for n in (1, 2, 3, 12)
            cmd = `$(MPI.mpiexec()) -n $n $(Base.julia_cmd()) $(joinpath(@__DIR__, "test_mpi.jl"))`
            p = run(ignorestatus(cmd))
            @test success(p)
        end
    end

    @testset "faces" begin
        dims = (16, 16, 16)
        core = CartesianIndices(map(n -> 2:n, dims))
        domain = CartesianIndices(map(n -> 1:(n + 1), dims))
        sort(setdiff(domain, core)) == sort(vcat(reshape.(Topologies.faces(domain, core), Ref(:))...))

        dims = (17, 18, 19)
        core = CartesianIndices(map(n -> 4:(n - 1), dims))
        domain = CartesianIndices(map(n -> -2:(n + 3), dims))
        sort(setdiff(domain, core)) == sort(vcat(reshape.(Topologies.faces(domain, core), Ref(:))...))
    end
end

@test_nowarn MPI.Finalize()
