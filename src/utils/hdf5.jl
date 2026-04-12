using HDF5: File, Group

const WRITABLE = Union{AbstractFloat,Integer,AbstractString,Bool}

function deepwrite(::Union{File,Group}, path::AbstractString, x)
    error("cannot write value of type $(typeof(x)) at path $path")
end
deepwrite(::Union{File,Group}, ::AbstractString, ::Nothing) = nothing
deepwrite(h::Union{File,Group}, path::AbstractString, x::WRITABLE) = write(h, path, x)
function deepwrite(h::Union{File,Group}, path::AbstractString, x::Array{<:WRITABLE})
    !isempty(x) && write(h, path, x)
end
function deepwrite(
    h::Union{File,Group}, path::AbstractString, dict::Union{Tuple,NamedTuple,AbstractDict,AbstractVector}
)
    for (key, value) in pairs(dict)
        deepwrite(h, joinpath(path, string(key)), value)
    end
end

deepread(h::Union{File,Group}, path::AbstractString) = deepread(h, path, h[path])
deepread(::Union{File,Group}, ::AbstractString, x) = read(x)
function deepread(h::Union{File,Group}, path::AbstractString, x::Union{File,Group})
    d = Dict{String,Any}()
    for key in keys(x)
        p = joinpath(path, key)
        d[key] = deepread(h, p, h[p])
    end
    d
end
