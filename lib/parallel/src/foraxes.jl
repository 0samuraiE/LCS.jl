@inline use_i32_if_possible(v::Integer) = v <= typemax(Int32) ? Int32(v) : v

@inline unwrap(x) = x
@inline unwrap(::Val{x}) where {x} = x

function _foraxes(f, ::KA.CPU, axes::Tuple{Vararg{AbstractUnitRange}})
    _foraxes(f, KA.CPU(), axes, Val(Threads.nthreads() == 1))
end

function _foraxes(f, ::KA.CPU, axes::Tuple{Vararg{AbstractUnitRange}}, ::Val{true})
    for i in CartesianIndices(axes)
        @inline f(Tuple(i)...)
    end
end

function _foraxes(f, ::KA.CPU, axes::Tuple{Vararg{AbstractUnitRange}}, ::Val{false})
    axes..., axis = axes
    Threads.@threads for j in axis
        for i in CartesianIndices(axes)
            @inline f(Tuple(i)..., j)
        end
    end
end

function _foraxes(
    f, backend::KA.GPU, axes::Union{NTuple{1,AbstractUnitRange},NTuple{2,AbstractUnitRange},NTuple{3,AbstractUnitRange}}
)
    i0 = use_i32_if_possible.(first.(axes))
    i1 = use_i32_if_possible.(last.(axes))
    (@kernel unsafe_indices = true function k()
        i = @index(Global, NTuple) .+ (i0 .- one(i0[1]))
        if all(i .<= i1)
            @inline f(i...)
        end
    end)(backend)(;
        ndrange=length.(axes)
    )
end

function _foraxes(f, backend::KA.GPU, axes::Tuple{Vararg{AbstractUnitRange}})
    caxes = CartesianIndices(axes)
    len = use_i32_if_possible(length(caxes))
    (@kernel unsafe_indices = true function k()
        i = @index(Global, Linear)
        if i <= len
            @inline f(Tuple(caxes[i])...)
        end
    end)(backend)(; ndrange=len)
end

"""
Abstract type for coloring schemes.
"""
abstract type Coloring end

"""
No coloring applied.
"""
struct NoColoring <: Coloring end

"""
    foraxes(f, backend, axes; coloring=NoColoring())

Dispatch parallel N-dimensional iteration over `axes` to `backend`.
"""
function foraxes(f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}; coloring::Coloring=NoColoring())
    foraxes_colored(f, backend, axes, coloring)
end

function foraxes_colored(f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}, ::NoColoring)
    _foraxes(f, backend, axes)
end
function foraxes_colored(f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}, coloring::Coloring)
    for color in (Val(0x00), Val(0x01))
        _foraxes_colored(f, backend, axes, coloring, color)
    end
end

"""
Red-black coloring; iterates cells where `sum(I) % 2 == color`.
"""
struct RedBlack <: Coloring end

@inline function _foraxes_colored(f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}, ::RedBlack, color)
    _foraxes(backend, axes) do i...
        sum(i) % 0x02 == unwrap(color) && @inline f(i...)
    end
end

"""
Optimized red-black coloring; halves iteration range in first dimension.
"""
struct RedBlackFast <: Coloring end

@inline function _foraxes_colored(f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}, ::RedBlackFast, color)
    a = use_i32_if_possible(first(axes[1]))
    b = use_i32_if_possible(last(axes[1]))
    n = use_i32_if_possible(length(axes[1]))
    halfaxes = (a:(a + cld(n, 2) - 1), Base.tail(axes)...)
    _foraxes(backend, halfaxes) do i, j...
        i2 = 0x02 * i - a
        i3 = i2 + oftype(i2, (sum((j...,); init=zero(i2)) - oftype(i2, a) + oftype(i2, unwrap(color))) & one(i2))
        if i3 <= oftype(i2, b)
            @inline f(i3, j...)
        end
    end
end

"""
    RedBlackBlock(size)

Block-based red-black coloring with block `size`.
"""
struct RedBlackBlock <: Coloring
    size :: Int
end

@inline function _foraxes_colored(
    f, backend::KA.Backend, axes::Tuple{Vararg{AbstractUnitRange}}, coloring::RedBlackBlock, color
)
    a = use_i32_if_possible.(first.(axes))
    b = use_i32_if_possible.(last.(axes))
    n = use_i32_if_possible.(length.(axes))
    bsize = use_i32_if_possible(coloring.size)
    block_axes = map(ni -> 1:cld(ni, bsize), n)
    _foraxes(backend, block_axes) do bi...
        if sum(bi) % 0x02 == unwrap(color)
            for li in CartesianIndices(map(_ -> 1:bsize, bi))
                is = map((ai, bii, lii) -> ai + bsize * (bii - 0x01) + lii - 0x01, a, bi, li.I)
                if all(map(<=, is, b))
                    @inline f(is...)
                end
            end
        end
    end
end
