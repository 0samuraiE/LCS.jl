struct Hasher{N,D<:Divisor}
    axes     :: NTuple{N,UnitRange{Int}}
    divisors :: NTuple{N,D}

    function Hasher(axes::NTuple{N,UnitRange{Int}}) where {N}
        lengths = map(length, axes)
        prod(lengths) <= typemax(UInt32) ||
            throw(ArgumentError("total elements $(prod(lengths)) exceeds typemax(UInt32)"))

        divisors = ntuple(i -> Divisor(lengths[i]), N)
        new{N,eltype(divisors)}(axes, divisors)
    end
end

@inline function encode(h::Hasher{N}, is::NTuple{N,Int}) where {N}
    sub2ind(h.axes, is)
end

@inline function decode(h::Hasher, hash::Int)
    ind2sub(h.axes, h.divisors, hash)
end

@inline function eachhash(h::Hasher)
    1:prod(length.(h.axes))
end

@inline @generated function sub2ind(ranges::NTuple{N,AbstractUnitRange}, I::NTuple{N,Int}) where {N}
    quote
        idx = 1
        stride = 1
        Base.@nexprs $N i -> begin
            offset = I[i] - first(ranges[i])
            idx += offset * stride
            stride *= length(ranges[i])
        end
        idx
    end
end

@inline @generated function ind2sub(
    ranges::NTuple{N,AbstractUnitRange}, divisors::NTuple{N,Divisor}, idx::Int
) where {N}
    quote
        offsets = map(first, ranges)
        linear = idx - 1

        Base.Cartesian.@ntuple $N i -> begin
            assume(0 <= linear <= typemax(Int32))

            if i == 1
                mod(linear, divisors[1]) + offsets[1]
            else
                prev = divisors[i - 1]
                linear = div(linear, prev)
                mod(linear, divisors[i]) + offsets[i]
            end
        end
    end
end
