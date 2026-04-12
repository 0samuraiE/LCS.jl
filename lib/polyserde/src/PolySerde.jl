module PolySerde
export @composite, @variant, @kind

const SerdeDict{V} = AbstractDict{String,V}
const SerdeDictType = Type{<:SerdeDict}

@nospecialize
struct Kind{T}
    Kind(T::Symbol) = new{T}()
end

struct Variant end
struct Composite end
struct Leaf end

compositetype(::Type) = Leaf()

"""
    kind2type(SuperType, Kind)

Maps kind symbol to concrete type. Registered by `@kind`.
"""
function kind2type end

"""
    type2kindstr(SuperType, instance)

Returns kind string for type instance.
"""
function type2kindstr end

"""
    deserialize(T, x; dicttype)

Converts dict or value to type `T`. Leaf types pass through, composite types
build from fields, variant types dispatch on `"kind"` field.

# Examples
```julia
# Leaf types (default)
deserialize(Int, 42)  # 42

# Composite types
struct Point
    x :: Int
    y :: Int
end
@composite Point
deserialize(Point, Dict("x" => 1, "y" => 2))  # Point(1, 2)

# Variant types
abstract type Shape end
@variant Shape

struct Circle <: Shape
    radius :: Float64
end
@composite Circle
@kind Shape "circle" Circle

deserialize(Shape, Dict("kind" => "circle", "params" => Dict("radius" => 1.0)))
# Circle(1.0)

# Passthrough for Tuple conversion
struct Grid
    dims :: Tuple{Int,Int,Int}
end
@composite Grid

function PolySerde.deserialize(::Composite, ::Type{Grid}, x::SerdeDict)
    passthrough = Dict("dims" => Tuple(x["dims"]))
    PolySerde.deserialize(Composite(), Grid, x, passthrough)
end

deserialize(Grid, Dict("dims" => [32, 64, 128]))  # Grid((32, 64, 128))
```

See also: [`serialize`](@ref), [`normalize`](@ref)
"""
deserialize(T::Type, x) = deserialize(compositetype(T), T, x)
deserialize(::Leaf, T::Type, x) = Base.convert(T, x)
deserialize(::Composite, T::Type, x::SerdeDict) = deserialize(Composite(), T, x, nothing)
function deserialize(::Composite, T::Type, x::SerdeDict, passthrough::Union{SerdeDict,Nothing})
    names = string.(fieldnames(T))

    miss = if isnothing(passthrough)
        setdiff(names, keys(x))
    else
        setdiff(names, union(keys(x), keys(passthrough)))
    end
    isempty(miss) || throw(ArgumentError("cannot deserialize $T, missing fields $(collect(miss))"))

    args = map(names, fieldtypes(T)) do name, type
        if !isnothing(passthrough) && haskey(passthrough, name)
            passthrough[name]
        else
            deserialize(type, x[name])
        end
    end
    T(args...)
end
function deserialize(::Variant, S::Type, x::SerdeDict)
    haskey(x, "kind") ||
        throw(ArgumentError("cannot deserialize to $S, expected key \"kind\", got $(collect(keys(x)))"))

    params = get(x, "params", Dict{String,Any}())
    deserialize(kind2type(S, Kind(Symbol(x["kind"]))), params)
end

"""
    serialize(T, x; dicttype)

Converts value to dict. Leaf types pass through, composite types map fields to
keys, variant types add `"kind"` and `"params"` fields.

# Examples
```julia
# Leaf types (default)
serialize(Int, 42)  # 42

# Composite types
struct Point
    x :: Int
    y :: Int
end
@composite Point
serialize(Point, Point(1, 2))  # Dict("x" => 1, "y" => 2)

# Variant types
abstract type Shape end
@variant Shape

struct Circle <: Shape
    radius :: Float64
end
@composite Circle
@kind Shape "circle" Circle

serialize(Shape, Circle(1.0))
# Dict("kind" => "circle", "params" => Dict("radius" => 1.0))
```

See also: [`deserialize`](@ref), [`@composite`](@ref), [`@variant`](@ref)
"""
serialize(T::Type, x; dicttype::SerdeDictType=Dict{String,Any}) = serialize(T, x, dicttype)
serialize(T::Type, x, dicttype::SerdeDictType) = serialize(compositetype(T), T, x, dicttype)
serialize(::Leaf, ::Type, x, dicttype::SerdeDictType) = x
serialize(::Composite, T::Type, x, dicttype::SerdeDictType) = serialize(Composite(), T, x, dicttype, nothing)
function serialize(::Composite, T::Type, x, dicttype::SerdeDictType, passthrough::Union{SerdeDict,Nothing})
    d = dicttype()
    for (name, type) in zip(fieldnames(T), fieldtypes(T))
        name_str = string(name)
        d[name_str] = if !isnothing(passthrough) && haskey(passthrough, name_str)
            passthrough[name_str]
        else
            serialize(type, getproperty(x, name), dicttype)
        end
    end
    d
end
function serialize(::Variant, S::Type, x, dicttype::SerdeDictType)
    d = dicttype()
    kindstr = type2kindstr(S, x)
    d["kind"] = kindstr
    params = serialize(kind2type(S, Kind(Symbol(kindstr))), x, dicttype)
    if !isempty(params)
        d["params"] = params
    end
    d
end

"""
    normalize(x; dicttype)

Recursively converts value to dict without type info or registration.

# Example
```julia
normalize((; x = 1, y = (; a = 2, b = 3)), Dict{String,Any})
# Dict("x" => 1, "y" => Dict("a" => 2, "b" => 3))
```

See also: [`deserialize`](@ref)
"""
normalize(x; dicttype::SerdeDictType=Dict{String,Any}) = normalize(x, dicttype)
normalize(x, dicttype::SerdeDictType) = normalize(compositetype(typeof(x)), x, dicttype)
normalize(::Leaf, x, dicttype::SerdeDictType) = x
normalize(::Leaf, x::Tuple, dicttype::SerdeDictType) = collect(map(y -> normalize(y, dicttype), x))
normalize(::Leaf, x::AbstractVector, dicttype::SerdeDictType) = map(x -> normalize(x, dicttype), x)
normalize(::Leaf, x::NamedTuple, dicttype::SerdeDictType) = normalize(Composite(), x, dicttype)
function normalize(::Leaf, d::AbstractDict, dicttype::SerdeDictType)
    dicttype(string(k) => normalize(v, dicttype) for (k, v) in d)
end
function normalize(::Composite, x, dicttype::SerdeDictType)
    d = dicttype()
    for name in fieldnames(typeof(x))
        name_str = string(name)
        d[name_str] = normalize(getproperty(x, name), dicttype)
    end
    d
end
function normalize(::Variant, x, dicttype::SerdeDictType)
    throw(ArgumentError("cannot normalize value of polymorphic supertype"))
end
@specialize

"""
    @composite Type

Enables serialization for struct with named fields.

See [`serialize`](@ref) and [`deserialize`](@ref) for usage examples.
"""
macro composite(type)
    ex = quote
        @nospecialize
        PolySerde.compositetype(::Type{<:$type}) = PolySerde.Composite()
        @specialize
    end
    esc(ex)
end

"""
    @variant SuperType
    @variant Union{Type1,Type2,...}

Enables polymorphic serialization using `"kind"` field to dispatch. Pair with
`@kind` and `@composite` on subtypes. Union form allows shared variants.

See [`serialize`](@ref) and [`deserialize`](@ref) for usage examples.
"""
macro variant(supertype)
    ex = quote
        @nospecialize
        PolySerde.compositetype(::Type{$supertype}) = PolySerde.Variant()
        @specialize
    end
    esc(ex)
end

"""
    @kind SuperType "discriminator" ConcreteType

Registers concrete type as variant. Discriminator appears in `"kind"` field.

See [`serialize`](@ref) and [`deserialize`](@ref) for usage examples.
"""
macro kind(supertype, name, subtype)
    ex = quote
        @nospecialize
        function PolySerde.kind2type(::Type{$supertype}, ::PolySerde.Kind{$(QuoteNode(Symbol(name)))})
            $subtype
        end
        PolySerde.type2kindstr(::Type{$supertype}, ::$subtype) = $name
        @specialize
    end
    esc(ex)
end
end
