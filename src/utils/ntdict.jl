const NTDict = Union{NamedTuple,Dict}

dict2nt(x) = x
dict2nt(x::Vector) = Tuple(dict2nt.(x))
dict2nt(x::AbstractDict) = NamedTuple((k => dict2nt(v) for (k, v) in x))

nt2dict(x; dicttype::Type{<:AbstractDict}=Dict) = nt2dict(x, dicttype)
nt2dict(x, ::Type{<:AbstractDict}) = x
nt2dict(x::Tuple, dicttype::Type{<:AbstractDict}) = [nt2dict(y, dicttype) for y in x]
nt2dict(x::NamedTuple, dicttype::Type{<:AbstractDict}) = dicttype((k => nt2dict(v, dicttype) for (k, v) in pairs(x)))
