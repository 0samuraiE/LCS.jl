module Offsets
export @offsets, @offsetviews, offset
using MacroTools

"""
    @offsets o ex

Offset array indices in `ex` by `o`. Use `\$A[i]` to escape. Supports integers, ranges,
`CartesianIndex`, `CartesianIndices`, and splatted tuples.

# Examples
```julia
@offsets 1 begin
    A[1, 1]                          # => A[2, 2]
    B[1:2, 3:4]                      # => B[2:3, 4:5]
    C[CartesianIndex(1, 2)]          # => C[CartesianIndex(2, 3)]
    D[CartesianIndices((1:2, 3:4))]  # => D[CartesianIndices((2:3, 4:5))]
    E[(1, 2:3)...]                   # => E[(2, 3:4)...]
    \$F[i, j]                        # => F[i, j] (escaped)
end
```
"""
macro offsets(o, ex)
    ret = MacroTools.prewalk(ex) do x
        if @capture(x, A_[is__])
            if Meta.isexpr(A, :$)
                :($(A.args[1])[$(is...)])
            else
                js = map(is) do i
                    if @capture(i, js_...)
                        :(Offsets.offset.($o, $js)...)
                    else
                        :(Offsets.offset($o, $i))
                    end
                end
                :($A[$(js...)])
            end
        else
            x
        end
    end
    esc(ret)
end

"""
    @offsetviews o ex

Equivalent to `@offsets o @views ex`.

# Examples
```julia
@offsetviews 5 begin
    A[i]    # => @views A[i + 5]
    B[i,j]  # => @views B[i + 5, j + 5]
end
```
"""
macro offsetviews(o, ex)
    esc(:(@offsets $o @views $ex))
end

"""
    offset(o, i::Integer)
    offset(o, i::CartesianIndex)
    offset(o, i::AbstractRange)
    offset(o, i::CartesianIndices)

Add offset `o` to index `i`.
"""
@inline offset(o::Integer, i::Integer) = i + o
@inline offset(o::Integer, i::CartesianIndex{N}) where {N} = i + CartesianIndex{N}(o)
@inline offset(o::Integer, i::AbstractRange) = i .+ o
@inline offset(o::Integer, i::CartesianIndices) = CartesianIndices(offset.(o, i.indices))
end # module Offsets
