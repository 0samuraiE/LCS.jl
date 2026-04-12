import AcceleratedKernels as AK

function sortperm!(perm::AbstractArray, A::AbstractArray, n::Integer=length(perm); temp=nothing)
    @views AK.sortperm!(perm[1:n], A[1:n]; temp=temp[1:n])
end

function sortbyperm!(A::AbstractArray, perms::AbstractArray, copy::AbstractArray)
    backend = KA.get_backend(A)
    Parallel.foraxes(backend, (eachindex(A, perms, copy),)) do i
        @inbounds begin
            copy[i] = A[perms[i]]
        end
    end
    # DO NOT REMOVE: synchronize() is required to prevent race conditions
    # between the parallel write operations and copyto!()
    KA.synchronize(backend)
    Base.copyto!(A, copy)
end

function sortbyperm!(A::AbstractArray, perms::AbstractArray, copy::AbstractArray, n::Integer)
    @views sortbyperm!(A[1:n], perms[1:n], copy[1:n])
end

function cumsum!(scan::AbstractArray, mask::AbstractArray, n::Integer=length(scan))
    @views Base.cumsum!(scan[1:n], mask[1:n])
end
