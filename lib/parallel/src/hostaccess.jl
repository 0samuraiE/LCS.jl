"""
    hostaccess(A, r)

Access elements `A[r]` on host; copies from GPU if needed.
"""
Base.@propagate_inbounds function hostaccess(A::AbstractArray, r)
    hostaccess(KA.get_backend(A), A, r)
end

Base.@propagate_inbounds function hostaccess(A::AbstractArray, r::Integer)
    hostaccess(KA.get_backend(A), A, r)[1]
end

Base.@propagate_inbounds function hostaccess(A::AbstractArray, rs...)
    hostaccess(KA.get_backend(A), A, rs...)
end

Base.@propagate_inbounds function hostaccess(::KA.CPU, A::AbstractVector, r)
    view(A, r)
end

Base.@propagate_inbounds function hostaccess(::KA.Backend, A::AbstractVector{T}, r) where {T}
    buf = zeros(T, length(r))
    copyto!(buf, view(A, r))
end

Base.@propagate_inbounds function hostaccess(::KA.CPU, A::AbstractArray, rs...)
    view(A, rs...)
end

Base.@propagate_inbounds function hostaccess(::KA.Backend, A::AbstractArray{T}, rs...) where {T}
    sub = view(A, rs...)
    buf = zeros(T, size(sub))
    copyto!(buf, sub)
    buf
end
