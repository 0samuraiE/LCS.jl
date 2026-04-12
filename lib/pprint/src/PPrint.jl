module PPrint
using Printf

abstract type PrintStyle end
struct Tree <: PrintStyle end
struct Default <: PrintStyle end

PrintStyle(::Type) = Default()
PrintStyle(::Type{<:NamedTuple}) = Tree()
PrintStyle(::Type{<:Tuple}) = Tree()

istree(::T) where {T} = istree(PrintStyle(T))
istree(::Type{T}) where {T} = istree(PrintStyle(T))
istree(::Default) = false
istree(::Tree) = true

rename(x) = string(x)
rename(x::Symbol) = rename(Val(x))
rename(::Val{x}) where {x} = string(x)

"""
    pprint(io::IO, x)

Pretty-print object `x` to `io` with tree-like formatting for nested structures.
"""
pprint(x) = pprint(stdout, x)
pprint(io::IO, x::T) where {T} = _pprint(io, PrintStyle(T), x, 0)
pprint(io::IO, x::Real) =
    let fmt = get(io, :realfmt, nothing)
        isnothing(fmt) ? print(io, x) : Printf.format(io, fmt, x)
    end

pprintln(x) = pprintln(stdout, x)
function pprintln(io::IO, x)
    pprint(io, x)
    println(io)
end

_pprint(io::IO, ::Default, x, level::Integer) = show(io, x)
function _pprint(io::IO, ::Tree, x::T, level::Integer) where {T}
    keys = fieldnames(T)
    span = maximum(k -> length(rename(k)), keys; init=0)

    first = true
    for key in keys
        value = getfield(x, key)
        isnothing(value) && continue
        is_empty_nt_recursive(value) && continue

        if !first
            print(io, "\n")
        end
        first = false

        renamed = rename(key)
        indent = "  "^level
        print(io, indent, renamed)

        if istree(value)
            print(io, " :\n")
            _pprint(io, Tree(), value, level + 1)
        else
            print(io, " "^(span - length(renamed)), " : ")
            pprint(io, value)
        end
    end
end

function _pprint(io::IO, ::Default, x::Vector, level::Integer)
    head = get(io, :vechead, 4)
    tail = get(io, :vectail, 4)
    print(io, "  "^level, "[")
    if length(x) <= head + tail
        for (i, item) in enumerate(x)
            pprint(io, item)
            i < length(x) && print(io, ", ")
        end
    else
        for i in 1:head
            pprint(io, x[i])
            i < head && print(io, ", ")
        end
        print(io, " ... ")
        for i in (length(x) - tail + 1):length(x)
            pprint(io, x[i])
            i < length(x) && print(io, ", ")
        end
    end
    print(io, "]")
end

is_empty_nt_recursive(x) = false
is_empty_nt_recursive(x::Union{Tuple,NamedTuple}) = isempty(x) || all(is_empty_nt_recursive, x)
end # module PPrint
