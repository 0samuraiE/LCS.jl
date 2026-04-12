function deepset!(d::Dict{T,Any}, value, keys::T...) where {T}
    current = d
    for k in keys[1:(end - 1)]
        current = get!(current, k, Dict{T,Any}())
    end
    current[keys[end]] = value
end
