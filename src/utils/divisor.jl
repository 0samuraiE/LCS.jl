# Simplified UInt32-only divisor for GPU/SIMD efficiency
# Trimmed down from Base.MultiplicativeInverses (stdlib)

abstract type Divisor end

struct Pow2Divisor <: Divisor
    shift   :: UInt32
    divisor :: UInt32
end

function Pow2Divisor(d::Integer)
    d = UInt32(d)
    ispow2(d) || throw(ArgumentError("Pow2Divisor requires power of 2, got $d"))
    Pow2Divisor(UInt32(trailing_zeros(d)), d)
end

@inline function Base.div(n::UInt32, pd::Pow2Divisor)
    n >> pd.shift
end

@inline function Base.mod(n::UInt32, pd::Pow2Divisor)
    n & (pd.divisor - 0x00000001)
end

@inline function Base.divrem(n::UInt32, pd::Pow2Divisor)
    q = n >> pd.shift
    r = n & (pd.divisor - 0x00000001)
    q, r
end

struct MagicDivisor <: Divisor
    multiplier :: UInt32
    shift      :: UInt32
    divisor    :: UInt32
end

function MagicDivisor(d::Integer)
    d = UInt32(d)

    m = UInt64(0)
    s = Int(0)
    for s in 0:62
        m = div((UInt64(1) << (32 + s)), d) + 1
        m <= typemax(UInt32) && break
    end
    MagicDivisor(UInt32(m), UInt32(s), d)
end

@inline function Base.div(n::UInt32, md::MagicDivisor)
    UInt32((UInt64(n) * UInt64(md.multiplier)) >> (32 + md.shift))
end

@inline function Base.mod(n::UInt32, md::MagicDivisor)
    n - div(n, md) * md.divisor
end

@inline function Base.divrem(n::UInt32, md::MagicDivisor)
    q = div(n, md)
    r = n - q * md.divisor
    q, r
end

@inline Base.div(n::Integer, d::Divisor) = div(UInt32(n), d)
@inline Base.mod(n::Integer, d::Divisor) = mod(UInt32(n), d)
@inline Base.divrem(n::Integer, d::Divisor) = divrem(UInt32(n), d)

function Divisor(d::Integer)
    d = UInt32(d)
    d == 1 && return Pow2Divisor(UInt32(0), UInt32(1))
    ispow2(d) ? Pow2Divisor(d) : MagicDivisor(d)
end
