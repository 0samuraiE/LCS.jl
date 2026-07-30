module AMDGPUExt
using Parallel
using AMDGPU

@inline function get_global_id(::NTuple{1,Any})
    (workitemIdx().x + (workgroupIdx().x - 0x01) * workgroupDim().x,)
end

@inline function get_global_id(::NTuple{2,Any})
    (
        workitemIdx().x + (workgroupIdx().x - 0x01) * workgroupDim().x,
        workitemIdx().y + (workgroupIdx().y - 0x01) * workgroupDim().y,
    )
end

@inline function get_global_id(::NTuple{3,Any})
    (
        workitemIdx().x + (workgroupIdx().x - 0x01) * workgroupDim().x,
        workitemIdx().y + (workgroupIdx().y - 0x01) * workgroupDim().y,
        workitemIdx().z + (workgroupIdx().z - 0x01) * workgroupDim().z,
    )
end

function Parallel._foraxes(
    f, ::ROCBackend, axes::Union{NTuple{1,AbstractUnitRange},NTuple{2,AbstractUnitRange},NTuple{3,AbstractUnitRange}}
)
    prod(length.(axes)) == 0 && return nothing

    i0 = Parallel.use_i32_if_possible.(first.(axes))
    i1 = Parallel.use_i32_if_possible.(last.(axes))
    function k()
        i = get_global_id(i0) .+ (i0 .- oftype(i0[1], 1))
        all(i .<= i1) || return nothing
        @inline f(i...)
        return nothing
    end
    k_ = @roc launch = false k()
    config = AMDGPU.launch_configuration(k_)
    groupsize = ntuple(d -> d == 1 ? config.groupsize : one(config.groupsize), length(axes))
    gridsize = cld.(length.(axes), groupsize)
    k_(; groupsize, gridsize)
end
end
