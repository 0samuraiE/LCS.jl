using Test

using FFTW
using HDF5
using LCS.Utils

@testset "concretize" begin
    let
        Utils.@concretize struct Example
            a :: Integer
            b :: AbstractString
            c :: AbstractFloat
        end

        x = Example(1, "test", 3.14)
        @test typeof(x) == Example{Int,String,Float64}
        @test x.a == 1
        @test x.b == "test"
        @test x.c == 3.14

        abstract type AbstractExample end
        Utils.@concretize struct Example2 <: AbstractExample
            a :: Integer
            b :: AbstractString
            c :: AbstractFloat
        end
        y = Example2(1, "test", 3.14)
        @test typeof(y) == Example2{Int,String,Float64}
        @test y.a == 1
        @test y.b == "test"
        @test y.c == 3.14

        @test supertype(Example) == Any
        @test supertype(Example2) == AbstractExample

        Utils.@kwconcretize struct Example3
            a :: Integer
            b :: AbstractString
            c :: AbstractFloat
        end

        z1 = Example3(; a=1, b="test", c=3.14)
        @test typeof(z1) == Example3{Int,String,Float64}
        z2 = Example3(1, "test", 3.14)
        @test typeof(z2) == Example3{Int,String,Float64}

        Utils.@concretize struct Example4{S<:Integer,T<:AbstractString,U<:AbstractFloat}
            a :: Tuple{S,T,U}
        end

        w = Example4((1, "test", 3.14))
        @test typeof(w) == Example4{Int,String,Float64,Tuple{Int,String,Float64}}
    end
end

@testset "deepset" begin
    let
        d = Dict{Symbol,Any}()
        Utils.deepset!(d, "val", :a, :b, :c)
        @test d[:a][:b][:c] == "val"

        Utils.deepset!(d, 42, :a, :b, :c)
        @test d[:a][:b][:c] == 42

        Utils.deepset!(d, "hello", :a, :x)
        @test d[:a][:x] == "hello"

        Utils.deepset!(d, true, :new, :path, :flag)
        @test d[:new][:path][:flag] == true

        d2 = Dict{String,Any}()
        Utils.deepset!(d2, 1, "a", "b")
        @test d2["a"]["b"] == 1
    end
end

@testset "fast" begin
    @testset "searchedsort2" begin
        let
            samples = [
                ([], 5),
                ([1], 0),
                ([1], 1),
                ([1], 2),
                ([1, 2, 2, 3, 5], 2),
                ([1, 2, 2, 3, 5], 4),
                ([1, 2, 2, 3, 5], 5),
                ([1, 2, 2, 3, 5], 6),
                ([1, 1, 1, 1, 1], 0),
                ([1, 1, 1, 1, 1], 1),
                ([1, 1, 1, 1, 1], 2),
            ]

            for (v, x) in samples
                @test Utils.searchsortedfirst(v, x) == searchsortedfirst(v, x)
                @test Utils.searchsortedlast(v, x) == searchsortedlast(v, x)
            end
        end
    end
end

@testset "hasher" begin
    @testset "sub/ind" begin
        ranges = (1:3, 1:4)
        divisors = Utils.Divisor.(length.(ranges))

        Utils.sub2ind(ranges, (1, 1)) == 1
        Utils.sub2ind(ranges, (3, 1)) == 3
        Utils.sub2ind(ranges, (1, 2)) == 4
        Utils.sub2ind(ranges, (3, 4)) == 12

        Utils.ind2sub(ranges, divisors, 1) == (1, 1)
        Utils.ind2sub(ranges, divisors, 3) == (3, 1)
        Utils.ind2sub(ranges, divisors, 4) == (1, 2)
        Utils.ind2sub(ranges, divisors, 12) == (3, 4)
    end

    let
        hasher = Utils.Hasher((0:5, 0:5, 0:5))
        x = Iterators.product(0:5, 0:5, 0:5)
        y = LinearIndices((6, 6, 6))

        @test all(Utils.encode.(Ref(hasher), x) .== y)
        @test all(Utils.decode.(Ref(hasher), y) .== x)
        @test all(Utils.decode.(Ref(hasher), Utils.encode.(Ref(hasher), x)) .== x)
        @test all(Utils.encode.(Ref(hasher), Utils.decode.(Ref(hasher), y)) .== y)
    end
end

@testset "ntdict" begin
    let
        @test Utils.dict2nt("test") == "test"
        @test Utils.nt2dict("test") == "test"

        @test Utils.dict2nt([1, 2, 3]) == (1, 2, 3)
        @test Utils.nt2dict((1, 2, 3)) == [1, 2, 3]

        d = Dict(:a => 1, :b => "hello", :c => [2, 3])
        nt = Utils.dict2nt(d)
        @test nt.a == 1
        @test nt.b == "hello"
        @test nt.c == (2, 3)

        d2 = Utils.nt2dict(nt)
        @test d2 == d

        d_nested = Dict(:outer => Dict(:inner => "value"))
        nt_nested = Utils.dict2nt(d_nested)
        @test nt_nested.outer.inner == "value"
        d_nested2 = Utils.nt2dict(nt_nested)
        @test d_nested2 == d_nested
    end
end

@testset "hdf5" begin
    @testset "scalar types" begin
        let
            tmp = tempname()
            x = Dict("int" => 42, "float" => 3.14, "string" => "hello", "bool" => true)

            h5open(tmp, "w") do h
                Utils.deepwrite(h, "data", x)
            end

            y = h5open(tmp, "r") do h
                Utils.deepread(h, "data")
            end

            @test y["int"] == 42
            @test y["float"] == 3.14
            @test y["string"] == "hello"
            @test y["bool"] == true
        end
    end

    @testset "array types" begin
        let
            tmp = tempname()
            x = Dict(
                "int_array" => [1, 2, 3, 4],
                "float_array" => [1.0, 2.0, 3.0],
                "string_array" => ["foo", "bar", "baz"],
                "bool_array" => [true, false, true],
            )

            h5open(tmp, "w") do h
                Utils.deepwrite(h, "data", x)
            end

            y = h5open(tmp, "r") do h
                Utils.deepread(h, "data")
            end

            @test y["int_array"] == [1, 2, 3, 4]
            @test y["float_array"] == [1.0, 2.0, 3.0]
            @test y["string_array"] == ["foo", "bar", "baz"]
            @test y["bool_array"] == [true, false, true]
        end
    end

    @testset "nested structures" begin
        let
            tmp = tempname()
            x = Dict("a" => 1, "b" => Dict("c" => 2, "d" => 3), "e" => [4, 5, 6])

            h5open(tmp, "w") do h
                Utils.deepwrite(h, "data", x)
            end

            y = h5open(tmp, "r") do h
                Utils.deepread(h, "data")
            end

            @test y == x
        end
    end

    @testset "NamedTuple" begin
        let
            tmp = tempname()
            x = (a=1, b=(c=2, d=3), e=[4, 5, 6])

            h5open(tmp, "w") do h
                Utils.deepwrite(h, "data", x)
            end

            y = h5open(tmp, "r") do h
                Utils.deepread(h, "data")
            end

            @test y["a"] == x.a
            @test y["b"]["c"] == x.b.c
            @test y["b"]["d"] == x.b.d
            @test y["e"] == x.e
        end
    end

    @testset "mixed types comprehensive" begin
        let
            tmp = tempname()
            x = Dict(
                "scalars" => Dict(
                    "int8" => Int8(1),
                    "int16" => Int16(2),
                    "int32" => Int32(3),
                    "int64" => Int64(4),
                    "float32" => Float32(1.5),
                    "float64" => Float64(2.5),
                ),
                "arrays" =>
                    Dict("vec_int" => [1, 2, 3], "vec_float" => [1.1, 2.2, 3.3], "vec_string" => ["a", "b", "c"]),
                "nested" => Dict("level1" => Dict("level2" => Dict("value" => 42))),
            )

            h5open(tmp, "w") do h
                Utils.deepwrite(h, "data", x)
            end

            y = h5open(tmp, "r") do h
                Utils.deepread(h, "data")
            end

            @test y == x
        end
    end
end

@testset "divisor" begin
    let
        pd = Utils.Pow2Divisor(8)
        @test div(16, pd) == 2
        @test mod(16, pd) == 0
        @test divrem(16, pd) == (2, 0)

        md = Utils.MagicDivisor(10)
        @test div(25, md) == 2
        @test mod(25, md) == 5
        @test divrem(25, md) == (2, 5)

        d = Utils.Divisor(4)
        @test div(15, d) == 3
        @test mod(15, d) == 3
        @test divrem(15, d) == (3, 3)

        d = Utils.Divisor(9)
        @test div(20, d) == 2
        @test mod(20, d) == 2
        @test divrem(20, d) == (2, 2)
    end
end
