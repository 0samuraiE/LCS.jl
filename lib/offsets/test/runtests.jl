using Test

using Offsets

@testset "Offsets" begin
    @testset "offset" begin
        @test offset(1, 1) == 2
        @test offset(1, CartesianIndex(1, 1, 1)) == CartesianIndex(2, 2, 2)
        @test offset(1, 1:3) == 2:4
        @test offset(1, CartesianIndices((0:1, 0:1))) == CartesianIndices((1:2, 1:2))
    end

    @testset "@offsets" begin
        A = reshape(1:25, (5, 5))

        @offsets 1 begin
            @test $A[1, 1] == 1
            @test A[1, 1] == $A[2, 2]
            @test A[1:2, 3:4] == $A[2:3, 4:5]
            @test A[CartesianIndex(1, 2)] == $A[CartesianIndex(2, 3)]
            @test A[CartesianIndices((1:2, 3:4))] == $A[CartesianIndices((2:3, 4:5))]
            @test A[(1, 2:3)...] == $A[(2, 3:4)...]
        end

        B = rand(3, 3, 3, 3)
        o = 1
        indices = (1, CartesianIndex(1), 1:2, CartesianIndices((1:2,)))
        @offsets o begin
            B[indices...] == $B[offset.(o, indices)...]
        end
    end
end
