function write_field(h::Union{HDF5.File,HDF5.Group}, name::String, A::AbstractArray, topo::Topologies.Topology)
    cart_rank = Topologies.cart_rank(topo)
    proc_dims = Topologies.proc_dims(topo)

    dims_l = size(A)
    dims_g = dims_l .* proc_dims

    indices = local_field_indices(dims_l, cart_rank)

    data = _create_dataset(h, name, eltype(A), dims_g, topo)
    data[indices...] = A
end

function read_field(h::Union{HDF5.File,HDF5.Group}, name::String, A::AbstractArray, topo::Topologies.Topology)
    cart_rank = Topologies.cart_rank(topo)
    proc_dims = Topologies.proc_dims(topo)

    dims_l = size(A)
    dims_g = dims_l .* proc_dims

    indices = local_field_indices(dims_l, cart_rank)

    data = _get_dataset(h, name, topo)
    dims_g == size(data) || throw(ArgumentError("dataset $name size mismatch, expected $dims_g got $(size(data))"))

    copyto!(A, data[indices...])
end

function local_field_indices(local_dims::NTuple{N,Integer}, cart_rank::NTuple{N,Integer}) where {N}
    map(local_dims, cart_rank) do n, i
        (1:n) .+ i * n
    end
end

function write_property(h::Union{HDF5.File,HDF5.Group}, name::String, A::AbstractVector, topo::Topologies.Topology)
    length(A) == 0 && return nothing

    rank = Topologies.linear_rank(topo)

    counts = Topologies.allgather(length(A), topo)
    count_g = sum(counts)

    data = _create_dataset(h, name, eltype(A), (count_g,), topo)
    indices = local_property_indices(counts, rank)

    data[indices] = A
end

function read_property(h::Union{HDF5.File,HDF5.Group}, name::String, A::AbstractVector, topo::Topologies.Topology)
    length(A) == 0 && return nothing

    rank = Topologies.linear_rank(topo)

    counts = Topologies.allgather(length(A), topo)
    count_g = sum(counts)

    indices = local_property_indices(counts, rank)

    data = _get_dataset(h, name, topo)
    length(data) == count_g ||
        throw(ArgumentError("dataset $name size mismatch, expected $count_g got $(length(data))"))

    copyto!(A, data[indices])
end

function local_property_indices(counts::AbstractVector{Int}, rank::Integer)
    start = @views sum(counts[1:(rank)]) + 1
    stop = @views sum(counts[1:(rank + 1)])
    start:stop
end

function _h5open(f, filename::AbstractString, mode::AbstractString, topo::Topologies.Topology)
    if Topologies.is_multi_processing(topo)
        comm = Topologies.comm()
        info = MPI.Info()
        h5open(f, filename, mode, comm, info)
    else
        h5open(f, filename, mode)
    end
end

function _get_dataset(h::Union{HDF5.File,HDF5.Group}, path::AbstractString, topo::Topologies.Topology)
    if Topologies.is_multi_processing(topo)
        h[path, dxpl_mpio=:collective]
    else
        h[path]
    end
end

function _create_dataset(
    h::Union{HDF5.File,HDF5.Group},
    path::AbstractString,
    dtype::Type,
    dspace_dims::Tuple{Vararg{Integer}},
    topo::Topologies.Topology,
)
    if Topologies.is_multi_processing(topo)
        create_dataset(h, path, dtype, dspace_dims; dxpl_mpio=:collective)
    else
        create_dataset(h, path, dtype, dspace_dims)
    end
end
