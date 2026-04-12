module CUDAExt
using CUDA
using KernelAbstractions
using PolySerde

@kind KernelAbstractions.Backend "cuda" CUDABackend
PolySerde.deserialize(::Type{CUDABackend}, ::PolySerde.SerdeDict) = CUDABackend()
PolySerde.serialize(::Type{CUDABackend}, ::CUDABackend, dicttype::PolySerde.SerdeDictType) = dicttype()
end
