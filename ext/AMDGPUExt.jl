module AMDGPUExt
using AMDGPU
using KernelAbstractions
using LCS
using Parallel
using PolySerde

PolySerde.@kind KernelAbstractions.Backend "amdgpu" AMDGPU.ROCBackend
PolySerde.deserialize(::Type{AMDGPU.ROCBackend}, ::PolySerde.SerdeDict) = AMDGPU.ROCBackend()
function PolySerde.serialize(
    ::Type{AMDGPU.ROCBackend}, ::AMDGPU.ROCBackend, dicttype::PolySerde.SerdeDictType
)
    dicttype()
end

# MIOpen does not support Float64 pooling, so use a backend-independent kernel
# over the coarse cells instead of NNlib.meanpool! on ROCArray.
function LCS.Flows.boxmean_l!(
    backend::AMDGPU.ROCBackend,
    U_coarse_l::AbstractArray,
    U::AbstractArray,
    halo_size::Integer,
    filter_dims::NTuple{3,Integer},
)
    kx, ky, kz = filter_dims
    w = 1 / LCS.FP(kx * ky * kz)

    domain = map(n -> 1:n, size(U_coarse_l))
    Parallel.foraxes(backend, domain) do ic, jc, kc
        @inbounds begin
            s = zero(eltype(U))
            for dk in 1:kz, dj in 1:ky, di in 1:kx
                i = halo_size + (ic - 1) * kx + di
                j = halo_size + (jc - 1) * ky + dj
                k = halo_size + (kc - 1) * kz + dk
                s += U[i, j, k]
            end
            U_coarse_l[ic, jc, kc] = s * w
        end
    end
end
end
