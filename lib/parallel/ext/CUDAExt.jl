module CUDAExt
using Parallel
using CUDA

@inline function get_global_id(::NTuple{1,Any})
    (threadIdx().x + (blockIdx().x - 0x01) * blockDim().x,)
end

@inline function get_global_id(::NTuple{2,Any})
    (threadIdx().x + (blockIdx().x - 0x01) * blockDim().x, threadIdx().y + (blockIdx().y - 0x01) * blockDim().y)
end

@inline function get_global_id(::NTuple{3,Any})
    (
        threadIdx().x + (blockIdx().x - 0x01) * blockDim().x,
        threadIdx().y + (blockIdx().y - 0x01) * blockDim().y,
        threadIdx().z + (blockIdx().z - 0x01) * blockDim().z,
    )
end

function Parallel._foraxes(
    f, ::CUDABackend, axes::Union{NTuple{1,AbstractUnitRange},NTuple{2,AbstractUnitRange},NTuple{3,AbstractUnitRange}}
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
    k_ = @cuda always_inline = true launch = false k()
    config = launch_configuration(k_.fun)
    threads = ntuple(d -> d == 1 ? config.threads : 1, length(axes))
    blocks = cld.(length.(axes), threads)
    k_(; threads, blocks)
end
end
