module ReferenceTests
export @test_reference

using JLD2
using LinearAlgebra
using Statistics
using Test

const ENV_JULIA_REFERENCETESTS_UPDATE = "JULIA_REFERENCETESTS_UPDATE"

"""
    @test_reference(file, dict, by=(==))

Test values against stored references or update them if `JULIA_REFERENCETESTS_UPDATE=true`.

The comparison function `by` is called as `by(reference, actual)`.
"""
macro test_reference(file, dict, kwargs...)
    quote
        if Base.get_bool_env(ReferenceTests.ENV_JULIA_REFERENCETESTS_UPDATE, false)
            ReferenceTests.update_reference($(esc(file)), $(esc(dict)))
            @info string("Updated reference file: ", $(esc(file)), ", keys: ", keys($(esc(dict))))
        else
            @test ReferenceTests.test_reference($(esc(file)), $(esc(dict)); $(map(esc, kwargs)...))
        end
    end
end

function report_mismatch(k::AbstractString, ref::AbstractArray, actual::AbstractArray)
    if size(ref) != size(actual)
        @warn "mismatch for key \"$k\"" ref_size = size(ref) actual_size = size(actual)
    else
        d = ref .- actual
        @warn "mismatch for key \"$k\"" mean = mean(abs, d) findmax = findmax(abs, d)
    end
end

function report_mismatch(k::AbstractString, ref, actual)
    @warn "mismatch for key \"$k\"" ref actual
end

function test_reference(file::AbstractString, dict::AbstractDict; by=(==))
    jldopen(file, "r") do f
        ok = true
        for (k, v) in dict
            haskey(f, k) || throw(ArgumentError("key $k not found in reference file $file"))
            r = f[k]
            if !by(r, v)
                report_mismatch(k, r, v)
                ok = false
            end
        end
        ok
    end
end

function update_reference(file::AbstractString, dict::AbstractDict)
    jldopen(file, "a+") do f
        for (k, v) in dict
            haskey(f, k) && delete!(f, k)
            f[k] = v
        end
    end
    jldrepack(file)
end

function jldrepack(src::AbstractString)
    tmp = tempname()
    try
        jldopen(tmp, "w") do fdst
            jldopen(src, "r") do fsrc
                _slim(fdst, fsrc)
            end
        end
        mv(tmp, src; force=true)
    catch
        isfile(tmp) && rm(tmp)
        rethrow()
    end
end

function _slim(fdst::Union{JLD2.JLDFile,JLD2.Group}, fsrc::Union{JLD2.JLDFile,JLD2.Group})
    for k in keys(fsrc)
        v = fsrc[k]

        if v isa JLD2.Group
            gdst = JLD2.Group(fdst, k)
            _slim(gdst, v)
        else
            fdst[k] = v
        end
    end
end
end # module ReferenceTests
