module FourierTools
using StaticArrays

"""
    fftindex(k, n)

Map FFT wavenumber `k` to a 1-based index for length `n`.
"""
fftindex(k::Integer, n::Integer) = k >= 0 ? k + 1 : n + k + 1

"""
    fftcoeff(A_hat, kk; normalize=true)

Get FFT coefficient at wavevector `kk`.
"""
function fftcoeff(A_hat::AbstractArray{<:Complex,N}, kk::NTuple{N,<:Integer}; normalize=true) where {N}
    i = fftindex.(kk, size(A_hat))
    if normalize
        A_hat[i...] / prod(size(A_hat))
    else
        A_hat[i...]
    end
end

"""
    rfftcoeff(A_hat, kk; normalize=true)

Get real-FFT coefficient at `kk` (Hermitian symmetry for `kk[1] < 0`).
"""
function rfftcoeff(A_hat::AbstractArray{<:Complex,N}, kk::NTuple{N,<:Integer}; normalize=true) where {N}
    if kk[1] < 0
        conj(rfftcoeff(A_hat, map(-, kk)))
    else
        i = fftindex.(kk, size(A_hat))
        if normalize
            A_hat[i...] / prod((size(A_hat, 1) * 2 - 2, size(A_hat)[2:end]...))
        else
            A_hat[i...]
        end
    end
end

"""
    realprod(x, y)

Compute `real(x*y)`.
"""
@inline function realprod(x::Complex, y::Complex)
    xr, xi = real(x), imag(x)
    yr, yi = real(y), imag(y)
    xr * yr - xi * yi
end

"""
    kset(k0)

Generate canonical 3D wavevectors (half-space) in shell `k0 - 0.5 <= |k| < k0 + 0.5`.
"""
function kset(k0::Integer)
    K = Vector{NTuple{3,Int64}}()
    for kz in (-k0):k0, ky in (-k0):k0, kx in 0:k0
        kk = (kx, ky, kz)
        if (k0 - 0.5)^2 <= sum(abs2, kk) < (k0 + 0.5)^2
            c = canon(kk)
            if !(c in K)
                push!(K, c)
            end
        end
    end
    K
end

function canon(k::Tuple{Vararg{Integer}})
    # determine sign to make first non-zero component positive
    for a in k
        a > 0 && return k
        a < 0 && return map(-, k)
    end
    k
end

"""
    kcount(k0)

Count 3D wavevectors (full cube) in shell `k0 - 0.5 <= |k| < k0 + 0.5`.
"""
function kcount(k0::Integer)
    c = 0
    for kz in (-k0):k0, ky in (-k0):k0, kx in (-k0):k0
        r2 = kx^2 + ky^2 + kz^2
        if (k0 - 0.5)^2 <= r2 < (k0 + 0.5)^2
            c += 1
        end
    end
    c
end

"""
    enespe(U_hat, V_hat, W_hat, k0; hasnyquist=false)

Compute energy in shell `k0 - 0.5 <= |k| < k0 + 0.5`.
"""
function enespe(
    U_hat::AbstractArray{T,3}, V_hat::AbstractArray{T,3}, W_hat::AbstractArray{T,3}, k0::Integer; hasnyquist=false
) where {T<:Complex}
    e = zero(real(T))
    for kz in (-k0):k0, ky in (-k0):k0, kx in 0:k0
        r2 = kx^2 + ky^2 + kz^2
        if (k0 - 0.5)^2 <= r2 < (k0 + 0.5)^2
            w = kx == 0 || (hasnyquist && kx == size(U_hat, 1) - 1) ? 1 : 2
            e +=
                w * (
                    abs2(rfftcoeff(U_hat, (kx, ky, kz))) +
                    abs2(rfftcoeff(V_hat, (kx, ky, kz))) +
                    abs2(rfftcoeff(W_hat, (kx, ky, kz)))
                )
        end
    end
    e / 2
end

hasnyquist(n::Integer) = iseven(n)

"""
    isselfconj(kk, dims)

Check if `kk` is self-conjugate for dimensions `dims`.
"""
@inline function isselfconj(kk::NTuple{N,Integer}, dims::NTuple{N,Integer}) where {N}
    all(zip(kk, dims)) do (k, n)
        k == 0 || (hasnyquist(n) && abs(k) == div(n, 2))
    end
end

alloc_for_crop(::Type{T}, kmax::Integer) where {T<:Complex} = zeros(T, kmax + 1, 2 * kmax + 1, 2 * kmax + 1)

"""
    crop(A_hat, kmax; temp=nothing)

Utility to extract a cropped spectrum suitable for use with `sirdft`.
"""
function crop(
    A_hat::AbstractArray{T,3}, kmax::Integer; temp::Union{AbstractArray{T,3},Nothing}=nothing
) where {T<:Complex}
    Ac_hat = if isnothing(temp)
        alloc_for_crop(T, kmax)
    else
        size(temp) == (kmax + 1, 2 * kmax + 1, 2 * kmax + 1) ||
            throw(ArgumentError("temp must have size (kmax + 1, 2*kmax + 1, 2*kmax + 1), got $(size(temp))"))
        temp
    end
    for kz in (-kmax):kmax, ky in (-kmax):kmax, kx in 0:kmax
        kk = (kx, ky, kz)
        Ac_hat[fftindex.(kk, size(Ac_hat))...] = rfftcoeff(A_hat, kk; normalize=true)
    end
    SArray{Tuple{size(Ac_hat)...}}(Ac_hat)
end

function powi(x::Symbol, n::Integer)
    n == 0 && return 1
    n > 0 && return Expr(:call, :*, x, powi(x, n - 1))
    n < 0 && return Expr(:call, :conj, powi(x, -n))
end

Base.@propagate_inbounds @generated function sirdft(
    ::Val{N}, Ac_hat::SArray{S,T}, exp1i::T, exp1j::T, exp1k::T
) where {N,S,T<:Complex}
    pre = Expr(:block)
    K = kset(N)

    body = Expr(:block)
    push!(body.args, :(a = zero(real(T))))
    for (kx, ky, kz) in K
        x = powi(:exp1i, kx)
        y = powi(:exp1j, ky)
        z = powi(:exp1k, kz)

        # assumes no self-conjugate modes in actual use (e.g. (0,0,0) not included)
        push!(body.args, :(a += 2 * realprod(Ac_hat[$(fftindex.((kx, ky, kz), S.parameters)...)], $x * $y * $z)))
    end
    push!(body.args, :(return a))

    Expr(:block, pre, body)
end

"""
    sirdft(Ac_hat, exp1i, exp1j, exp1k, kmin, kmax)

Compute sparse inverse DFT over shells `kmin - 0.5 <= |k| < kmax + 0.5`.

Typically, `Ac_hat` is obtained by applying `crop` to a full spectrum `A_hat`.
"""
Base.@propagate_inbounds function sirdft(
    Ac_hat::SArray{S,T}, exp1i::T, exp1j::T, exp1k::T, kmin::Integer, kmax::Integer
) where {S,T}
    sirdft(Val(kmin), Val(kmax), Ac_hat, exp1i, exp1j, exp1k)
end
Base.@propagate_inbounds @generated function sirdft(
    ::Val{KMIN}, ::Val{KMAX}, Ac_hat::SArray{S,T}, exp1i::T, exp1j::T, exp1k::T
) where {KMIN,KMAX,S,T<:Complex}
    body = Expr(:block)
    push!(body.args, :(a = zero(real(T))))
    for K in KMIN:KMAX
        push!(body.args, :(a += sirdft(Val($K), Ac_hat, exp1i, exp1j, exp1k)))
    end
    push!(body.args, :(return a))
    body
end
end # module FourierTools
