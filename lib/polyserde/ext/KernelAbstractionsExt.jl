module KernelAbstractionsExt
using KernelAbstractions
using PolySerde

@variant KernelAbstractions.Backend
@kind KernelAbstractions.Backend "cpu" KernelAbstractions.CPU
PolySerde.deserialize(::Type{KernelAbstractions.CPU}, ::PolySerde.SerdeDict) = KernelAbstractions.CPU()
function PolySerde.serialize(
    ::Type{KernelAbstractions.CPU}, ::KernelAbstractions.CPU, dicttype::PolySerde.SerdeDictType
)
    dicttype()
end
end
