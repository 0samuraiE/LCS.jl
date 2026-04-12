module Schema
export generate

using MacroTools
using JSON
using OrderedCollections
using YAML

const DRAFT7_HEADER = OrderedDict("\$schema" => "http://json-schema.org/draft-07/schema#")
const DIR_META_SCHEMAS = joinpath(@__DIR__, "../schemas")
const FILES_META_SCHEMA_COMMON = "common.meta.yaml"
const FILES_META_SCHEMA_FLOW = "flow.meta.yaml"
const FILES_META_SCHEMA_PARTICLE = "particle.meta.yaml"
const FILE_SCHEMA = "lcs-schema.json"

"""
    generate()

Generate JSON Schema (Draft 7) for LCS configuration validation.

Produces `lcs-schema.json` that validates LCS configuration files. The schema
defines two simulation modes:
- `flow`: flow-only simulation
- `flow-particle`: combined flow and particle simulation

# Meta-schema notation

YAML meta-schemas provide a compact notation for defining JSON Schema
structures. The notation supports:

## Combinators
- `oneOf`, `allOf`, `anyOf`, `not`: schema composition
- Use `oneOf` with string constants instead of `enum`

## Types
- Primitive: `string`, `number`, `integer`, `boolean`
- Fixed-length arrays: `string(n)`, `number(n)`, `integer(n)`, `boolean(n)`

## Semantics
- Bare type strings declare types
- Non-type strings declare string constants
- Dictionaries declare objects (all fields required, no additional properties)
- Explicit `type: "array"` with `items` declares variable-length arrays
"""
function generate()
    dicttype = OrderedDict{String,Any}
    dict = dicttype()

    common = YAML.load_file(joinpath(DIR_META_SCHEMAS, FILES_META_SCHEMA_COMMON); dicttype)
    flow = YAML.load_file(joinpath(DIR_META_SCHEMAS, FILES_META_SCHEMA_FLOW); dicttype)
    particle = YAML.load_file(joinpath(DIR_META_SCHEMAS, FILES_META_SCHEMA_PARTICLE); dicttype)

    dict = dicttype(#
        "oneOf" => [
            dicttype("mode" => dicttype("kind" => "flow"), common..., "flow" => flow),
            dicttype(
                "mode" => dicttype("kind" => "flow-particle"),
                common...,
                "flow" => flow,
                "particles" => dicttype("type" => "array", "items" => particle),
            ),
        ],
    )

    schema = to_schema(dict)
    schema = dicttype(DRAFT7_HEADER..., schema...)
    write(FILE_SCHEMA, JSON.json(schema))
end

function match_type_string(s)
    m = match(r"(string|integer|number|boolean)(?:\((\d+)\))?", s)
    if m !== nothing
        type = m.captures[1]
        size = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
        (type, size)
    else
        nothing
    end
end

function to_schema(x::String)
    m = match_type_string(x)
    if isnothing(m)
        OrderedDict("type" => "string", "const" => x)
    else
        type, size = m
        if size == 1
            OrderedDict("type" => type)
        else
            OrderedDict(#
                "type" => "array",
                "items" => OrderedDict("type" => type),
                "minItems" => size,
                "maxItems" => size,
            )
        end
    end
end

function to_schema(x::Vector)
    map(to_schema, x)
end

function to_schema(x::OrderedDict)
    keys_ = keys(x)
    if length(keys_) == 1 && first(keys_) in ("oneOf", "allOf", "anyOf", "not")
        key = first(keys_)
        OrderedDict(key => to_schema(x[key]))
    elseif get(x, "type", "") == "array"
        OrderedDict(#
            x...,
            "items" => to_schema(x["items"]),
        )
    else
        OrderedDict(#
            "type" => "object",
            "properties" => OrderedDict(k => to_schema(v) for (k, v) in x),
            "required" => keys(x),
            "additionalProperties" => false,
        )
    end
end
end
