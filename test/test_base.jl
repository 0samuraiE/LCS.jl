using Test

using LCS
using OrderedCollections
using PolySerde
using Topologies
using YAML

@testset "mapping" begin
    let
        d1 = YAML.load_file(joinpath(@__DIR__, "../precompile.lcs-yaml"); dicttype=OrderedDict{String,Any})
        d1["resume"] = false
        config1 = PolySerde.deserialize(LCS.Config, d1)
        d2 = PolySerde.serialize(LCS.Config, config1; dicttype=OrderedDict{String,Any})
        @test d1 == d2
        config2 = PolySerde.deserialize(LCS.Config, d2)
        @test config1 == config2

        d3 = deepcopy(d1)
        d3["particles"] = nothing
        config3 = PolySerde.deserialize(LCS.Config, d3)
        d4 = PolySerde.serialize(LCS.Config, config3; dicttype=OrderedDict{String,Any})
        @test d3 == d4
        config4 = PolySerde.deserialize(LCS.Config, d4)
        @test config3 == config4
    end
end

@testset "config" begin
    let
        grid = LCS.Grid((16, 16, 16), 3)
        topo = Topologies.Mock.mocktopology(; linear_rank=3, proc_dims=(4, 2, 1))
        cart_rank = Topologies.cart_rank(topo)
        @test LCS.dims_g(grid) == (16, 16, 16)
        @test LCS.dims_l(grid, topo) == (4, 8, 16)
        @test LCS.spacings(grid) == LCS.DOMAIN_LENGTH ./ (16, 16, 16)
        @test LCS.lengths_g() == (LCS.DOMAIN_LENGTH, LCS.DOMAIN_LENGTH, LCS.DOMAIN_LENGTH)
        @test LCS.lengths_l(topo) == (LCS.DOMAIN_LENGTH / 4, LCS.DOMAIN_LENGTH / 2, LCS.DOMAIN_LENGTH)
        @test LCS.origins_g() == (LCS.DOMAIN_ORIGIN, LCS.DOMAIN_ORIGIN, LCS.DOMAIN_ORIGIN)
        @test LCS.origins_l(topo) == LCS.DOMAIN_LENGTH .* (1 / 4, 1 / 2, 0)
    end
end
