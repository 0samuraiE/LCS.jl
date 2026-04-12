using Test

using PolySerde
import PolySerde as PS

@testset "PolySerde" begin
    @testset "Kind" begin
        @test PS.Kind(:test) isa PS.Kind{:test}
    end

    @testset "serialize primitives" begin
        @test PS.serialize(Int, 5) == 5
        @test PS.serialize(String, "hello") == "hello"
        @test PS.serialize(Float64, 3.14) == 3.14
    end

    @testset "deserialize primitives" begin
        @test PS.deserialize(Int, 5) == 5
        @test PS.deserialize(String, "hello") == "hello"
        @test PS.deserialize(Vector, [1, 2, 3]) == [1, 2, 3]
    end

    @testset "@composite" begin
        struct TestStruct
            x :: Int
            y :: String
        end
        @composite TestStruct

        obj = TestStruct(10, "test")
        dict = Dict("x" => 10, "y" => "test")

        @test PS.serialize(TestStruct, obj) == dict
        @test PS.deserialize(TestStruct, dict) == obj
    end

    @testset "@variant" begin
        abstract type SuperType end
        @variant SuperType

        struct Type1 <: SuperType
            value :: Int
        end
        @composite Type1
        @kind SuperType "type1" Type1

        struct Type2 <: SuperType
            name  :: String
            count :: Int
        end
        @composite Type2
        @kind SuperType "type2" Type2

        obj1 = Type1(42)
        dict1 = Dict("kind" => "type1", "params" => Dict("value" => 42))

        @test PS.serialize(SuperType, obj1) == dict1
        @test PS.deserialize(SuperType, dict1) == obj1

        obj2 = Type2("example", 10)
        dict2 = Dict("kind" => "type2", "params" => Dict("name" => "example", "count" => 10))

        @test PS.serialize(SuperType, obj2) == dict2
        @test PS.deserialize(SuperType, dict2) == obj2
    end

    @testset "ignore extra fields" begin
        struct IgnoreTest
            x :: Int
        end
        @composite IgnoreTest

        dict = Dict("x" => 5, "y" => "should_be_ignored")
        obj = PS.deserialize(IgnoreTest, dict)
        @test obj.x == 5
    end

    @testset "missing required fields" begin
        struct RequiredTest
            x :: Int
            y :: String
        end
        @composite RequiredTest

        dict = Dict("x" => 5)
        @test_throws ArgumentError PS.deserialize(RequiredTest, dict)
    end

    @testset "leaf types round-trip" begin
        @test PS.deserialize(Int, PS.serialize(Int, 42)) == 42
        @test PS.deserialize(String, PS.serialize(String, "test")) == "test"
        @test PS.deserialize(Float64, PS.serialize(Float64, 3.14)) == 3.14
    end

    @testset "composite round-trip" begin
        struct Point
            x :: Int
            y :: Int
        end
        @composite Point

        obj = Point(1, 2)
        dict = Dict("x" => 1, "y" => 2)

        @test PS.deserialize(Point, PS.serialize(Point, obj)) == obj
        @test PS.serialize(Point, PS.deserialize(Point, dict)) == dict
    end

    @testset "variant round-trip" begin
        abstract type Shape end
        @variant Shape

        struct Circle <: Shape
            radius :: Float64
        end
        @composite Circle
        @kind Shape "circle" Circle

        obj = Circle(5.0)
        dict = Dict("kind" => "circle", "params" => Dict("radius" => 5.0))

        @test PS.deserialize(Shape, PS.serialize(Shape, obj)) == obj
        @test PS.serialize(Shape, PS.deserialize(Shape, dict)) == dict
    end

    @testset "nested composite round-trip" begin
        struct Inner
            value :: Int
        end
        @composite Inner

        struct Outer
            inner :: Inner
            name  :: String
        end
        @composite Outer

        obj = Outer(Inner(42), "test")
        dict = Dict("inner" => Dict("value" => 42), "name" => "test")

        @test PS.deserialize(Outer, PS.serialize(Outer, obj)) == obj
        @test PS.serialize(Outer, PS.deserialize(Outer, dict)) == dict
    end

    @testset "normalize" begin
        nt = (; x=1, y=(; a=2, b=3))
        @test PS.normalize(nt) == Dict("x" => 1, "y" => Dict("a" => 2, "b" => 3))

        @test PS.normalize((1, 2, 3)) == [1, 2, 3]
        @test PS.normalize([[1, 2], [3, 4]]) == [[1, 2], [3, 4]]

        mixed = (; coords=(1, 2, 3), meta=Dict(:name => "test"))
        normalized = PS.normalize(mixed)
        @test normalized["coords"] == [1, 2, 3]
        @test normalized["meta"] == Dict("name" => "test")
    end

    @testset "passthrough" begin
        struct Grid
            dims :: Tuple{Int,Int,Int}
        end
        @composite Grid

        function PS.deserialize(::PS.Composite, ::Type{Grid}, x::PS.SerdeDict)
            passthrough = Dict("dims" => Tuple(x["dims"]))
            PS.deserialize(PS.Composite(), Grid, x, passthrough)
        end

        dict = Dict("dims" => [32, 64, 128])
        obj = PS.deserialize(Grid, dict)
        @test obj.dims == (32, 64, 128)
    end

    @testset "variant empty params" begin
        abstract type Token end
        @variant Token

        struct EOF <: Token end
        @composite EOF
        @kind Token "eof" EOF

        result = PS.serialize(Token, EOF())
        @test result == Dict("kind" => "eof")
        @test !haskey(result, "params")

        @test PS.deserialize(Token, Dict("kind" => "eof")) == EOF()
    end
end
