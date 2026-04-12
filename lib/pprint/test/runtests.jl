using Test

using Printf
using PPrint: PPrint, PrintStyle, istree, Tree, Default, pprint

@testset "PPrint" begin
    @test PrintStyle(Int) == Default()
    @test PrintStyle(Float64) == Default()
    @test PrintStyle(NamedTuple) == Tree()

    @test istree(1) == false
    @test istree((; x=1)) == true

    @test istree(Int) == false
    @test istree(NamedTuple) == true

    struct Summary
        mean
        max
        min
    end

    function PPrint.pprint(io::IO, s::Summary)
        print(io, "(")
        pprint(io, s.mean)
        print(io, ", ")
        pprint(io, s.max)
        print(io, ", ")
        pprint(io, s.min)
        print(io, ")")
    end

    struct Inner2
        e
        f
    end

    struct Inner
        c
        d
        inner2 :: Inner2
    end
    PPrint.PrintStyle(::Type{<:Inner}) = Tree()

    struct Outer
        a
        b
        c     :: Summary
        inner :: Inner
    end
    PPrint.PrintStyle(::Type{<:Outer}) = Tree()

    obj = Outer(1.0, 2 + 3im, Summary(2, 3, 1), Inner(pi, 4.0, Inner2(true, 12)))

    io = IOBuffer()
    pprint(io, obj)
    ref = """
          a     : 1.0
          b     : 2 + 3im
          c     : (2, 3, 1)
          inner :
            c      : π
            d      : 4.0
            inner2 : Inner2(true, 12)"""
    @test String(take!(io)) == ref

    @testset "Vector" begin
        let v = [1.0, 2.0, 3.0]
            io = IOBuffer()
            pprint(io, v)
            @test String(take!(io)) == "[1.0, 2.0, 3.0]"
        end

        let v = collect(1.0:10.0)
            io = IOBuffer()
            pprint(io, v)
            @test String(take!(io)) == "[1.0, 2.0, 3.0, 4.0 ... 7.0, 8.0, 9.0, 10.0]"
        end

        let v = collect(1.0:8.0)
            io = IOBuffer()
            pprint(io, v)
            @test String(take!(io)) == "[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]"
        end

        let v = collect(1.0:10.0)
            io = IOContext(IOBuffer(), :vechead => 2, :vectail => 1)
            pprint(io, v)
            @test String(take!(io.io)) == "[1.0, 2.0 ... 10.0]"
        end
    end

    @testset "realfmt" begin
        let io = IOBuffer()
            pprint(io, 3.14)
            @test String(take!(io)) == "3.14"
        end

        let io = IOContext(IOBuffer(), :realfmt => Printf.Format("%+.2E"))
            pprint(io, 3.14)
            @test String(take!(io.io)) == "+3.14E+00"
        end
    end
end
