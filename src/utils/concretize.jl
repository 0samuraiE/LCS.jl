using MacroTools

function _build_struct(SName, SuperType, fields, wheres; kwdef)
    parsed = map(fields) do f
        @capture(f, name_::T_) || error("Field must have the form `name::Type`, got `$f`")
        (name, T)
    end

    ConcreteTypes = []
    for w in wheres
        if @capture(w, ConcreteType_ <: AbstractType_)
            push!(ConcreteTypes, ConcreteType)
        else
            push!(ConcreteTypes, w)
        end
    end

    fields = []
    for (name, AbstractType) in parsed
        if AbstractType in ConcreteTypes
            push!(fields, :($name::$AbstractType))
        else
            ConcreteType = Symbol(name, "_T")
            push!(wheres, :($ConcreteType <: $AbstractType))
            push!(fields, :($name::$ConcreteType))
        end
    end

    if kwdef
        quote
            Base.@__doc__ @kwdef struct $SName{$(unblock.(wheres)...)} <: $SuperType
                $(unblock.(fields)...)
            end
        end
    else
        quote
            Base.@__doc__ struct $SName{$(unblock.(wheres)...)} <: $SuperType
                $(unblock.(fields)...)
            end
        end
    end
end

function _build_show(SName)
    quote
        function Base.show(io::IO, ::MIME"text/plain", x::$SName)
            keys = fieldnames(typeof(x))
            values = getfield.(Ref(x), keys)
            max_length = maximum(length, String.(keys))
            println(io, "Concretized struct ", $SName, " <: ", supertype(typeof(x)))
            map(keys, values) do key, value
                println(io, "    ", rpad(key, max_length), " :: ", typeof(value))
            end
            println(io, "end")
        end
    end
end

function _parse_expr(ex)
    #! format: off
    @capture(
        ex,
        (struct SName_{wheres__} <: SuperType_
            fields__
        end) |
        (struct SName_ <: SuperType_
            fields__
        end) |
        (struct SName_{wheres__}
            fields__
        end) |
        (struct SName_
            fields__
        end),
        ) || error("Expected struct definition, got `$ex`")
    #! format: on
    if isnothing(SuperType)
        SuperType = Any
    end
    if isnothing(wheres)
        wheres = []
    end
    SName, SuperType, fields, wheres
end

macro concretize(ex)
    SName, SuperType, fields, wheres = _parse_expr(ex)
    Expr(:block, esc(_build_struct(SName, SuperType, fields, wheres; kwdef=false)), esc(_build_show(SName)))
end

macro kwconcretize(ex)
    SName, SuperType, fields, wheres = _parse_expr(ex)
    Expr(:block, esc(_build_struct(SName, SuperType, fields, wheres; kwdef=true)), esc(_build_show(SName)))
end
