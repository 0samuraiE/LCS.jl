# performant reshape
# https://github.com/JuliaLang/julia/issues/36313

function reshape(v, dims...)
    Base.reshape(view(v, :), dims...)
end

function searchsortedfirst(v::Union{AbstractVector{T},Tuple{Vararg{T}}}, x::T) where {T}
    lo = 1
    hi = length(v) + 1
    len = hi - lo
    while len != 0
        half = len >>> 1
        m = lo + half
        isless = v[m] < x
        lo = isless * (m + 1) + (!isless) * lo
        len = isless * (len - (half + 1)) + (!isless) * half
    end
    lo
end

function searchsortedfirst(v::Union{AbstractVector{T},Tuple{Vararg{T}}}, x::T, n::Integer) where {T}
    searchsortedfirst(view(v, 1:n), x)
end

function searchsortedlast(v::Union{AbstractVector{T},Tuple{Vararg{T}}}, x::T) where {T}
    lo = 0
    hi = length(v) + 1
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 1)
        isless = x < v[m]
        hi = isless * m + (!isless) * hi
        lo = (!isless) * m + isless * lo
    end
    lo
end

function searchsortedlast(v::Union{AbstractVector{T},Tuple{Vararg{T}}}, x::T, n::Integer) where {T}
    searchsortedlast(view(v, 1:n), x)
end

function norm(x)
    s = sum(_x -> _x^2, x)
    assume(s >= 0)
    sqrt(s)
end

# assume typemin(T) <= x <= typemax(T)
function unsafe_floor(::Type{T}, x::Real) where {T<:Integer}
    unsafe_trunc(T, Base.floor(x))
end
unsafe_floor(x::Real) = unsafe_floor(Int32, x)

# Wrap x into [O, O+L) using floor; faster and branch-free vs mod-based form.
function wrap(x::Real, O::Real, L::Real)
    x - L * Base.floor((x - O) / L)
end

function sliceshape(A::AbstractArray, dims::Integer...)
    Utils.reshape(view(A, 1:prod(dims)), dims...)
end
