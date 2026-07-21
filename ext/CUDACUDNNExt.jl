module CUDACUDNNExt
using CUDA
using cuDNN
using LCS
using NNlib

function LCS.Flows.boxmean_l!(
    ::CUDA.CUDABackend,
    U_coarse_l::AbstractArray,
    U::AbstractArray,
    halo_size::Integer,
    filter_dims::NTuple{3,Integer},
)
    U = LCS.Utils.reshape(U, size(U)..., 1, 1)
    U_coarse_l = LCS.Utils.reshape(U_coarse_l, size(U_coarse_l)..., 1, 1)

    k = stride = filter_dims
    padding = -halo_size
    meanpool!(U_coarse_l, U, PoolDims(size(U), k; stride, padding))
end
end
