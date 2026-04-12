const N_COMMS = 2
const OFFSET = 1000

"""
    HaloBuffer

Abstract type for halo exchange buffers.
"""
abstract type HaloBuffer end

"""
    SPHaloBuffer <: HaloBuffer

Single-process halo buffer for non-distributed execution.
"""
struct SPHaloBuffer <: HaloBuffer end

"""
    MPHaloBuffer{V<:AbstractVector{<:Real}} <: HaloBuffer

Multi-process halo buffer for distributed execution.
"""
struct MPHaloBuffer{V<:AbstractVector{<:Real}} <: HaloBuffer
    requests :: MPI.MultiRequest
    sendbufs :: NTuple{Topologies.N_DIRECTIONS,V}
    recvbufs :: NTuple{Topologies.N_DIRECTIONS,V}
    nwrites  :: Vector{Int}
end

function HaloBuffer(
    backend::KA.Backend, topo::Topology, ::Type{T}, dims::NTuple{N_DIMS,Integer}, halo_size::Integer, narrays::Integer
) where {T<:Real}
    requests = MPI.MultiRequest(Topologies.N_DIRECTIONS * N_COMMS)
    sendbufs = _allocate_comm_bufs(backend, T, dims, halo_size, narrays)
    recvbufs = _allocate_comm_bufs(backend, T, dims, halo_size, narrays)
    nwrites = zeros(Int, Topologies.N_DIRECTIONS)

    if is_multi_processing(topo)
        MPHaloBuffer(requests, sendbufs, recvbufs, nwrites)
    else
        SPHaloBuffer()
    end
end

function _allocate_comm_bufs(
    backend::KA.Backend, ::Type{T}, dims::NTuple{N_DIMS,Integer}, halo_size::Integer, narrays::Integer
) where {T<:Real}
    map(Topologies.DIRECTIONS) do direction
        idst = recv_indices(direction, dims, halo_size)
        len = narrays * length(CartesianIndices(idst))
        KA.zeros(backend, T, len)
    end
end

"""
    HaloStyle

Abstract type for halo exchange styles.
"""
abstract type HaloStyle end

"""
    FaceHalo <: HaloStyle

Face-only halo exchange style (6 directions).

# Fields
- `size::Int`: Halo region size
"""
struct FaceHalo <: HaloStyle
    size :: Int
end
comm_directions(::FaceHalo) = Topologies.DIRECTIONS[1:6]
halo_size_comm(h::FaceHalo) = h.size

"""
    FullHalo <: HaloStyle

Full halo exchange style (26 directions).

# Fields
- `size::Int`: Halo region size
"""
struct FullHalo <: HaloStyle
    size :: Int
end
comm_directions(::FullHalo) = Topologies.DIRECTIONS
halo_size_comm(h::FullHalo) = h.size

"""
    synchalo!(style, As, buf, dims, halo_size, backend, topo)

Synchronize halo regions across all processes.

Perform halo exchange for arrays `As` using the specified style and buffer.
The operation is synchronized across all processes in the topology.
"""
function synchalo!(
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    buf::HaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    synchalo!(() -> nothing, style, As, buf, dims, halo_size, backend, topo)
end
function synchalo!(
    compute::Function,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    buf::HaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    packhalo!(processing(topo), style, As, buf, dims, halo_size, backend, topo)
    commhalo!(processing(topo), style, As, buf, dims, halo_size, backend, topo)
    t = @async waithalo(processing(topo), buf)
    compute()
    wait(t)
    unpackhalo!(processing(topo), style, As, buf, dims, halo_size, backend, topo)
end

function packhalo!(
    ::SingleProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    ::SPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
) end

function commhalo!(
    ::SingleProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    ::SPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    halo_size_comm = Topologies.halo_size_comm(style)
    directions = comm_directions(style)

    for A in As, direction in directions
        irecv = offset_recv_indices(direction, dims, halo_size_comm, halo_size)
        isend = offset_send_indices_to_self(direction, dims, halo_size_comm, halo_size)
        @views recv, send = A[irecv...], A[isend...]
        copyto!(recv, send)
    end

    KA.synchronize(backend)
end

function waithalo(::SingleProcessing, ::SPHaloBuffer) end

function unpackhalo!(
    ::SingleProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    ::SPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
) end

function packhalo!(
    ::MultiProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    buf::MPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    (; sendbufs, nwrites) = buf

    halo_size_comm = Topologies.halo_size_comm(style)
    directions = comm_directions(style)
    neighbors = topo.neighbor_linear_ranks
    self = Topologies.linear_rank(topo)

    for (i, (direction, neighbor, sendbuf)) in enumerate(zip(directions, neighbors, sendbufs))
        nwrite = 0
        for A in As
            if neighbor == self
                irecv = offset_recv_indices(direction, dims, halo_size_comm, halo_size)
                isend = offset_send_indices_to_self(direction, dims, halo_size_comm, halo_size)
                @views recv, send = A[irecv...], A[isend...]
                copyto!(recv, send)
            else
                isend = offset_send_indices_to_other(direction, dims, halo_size_comm, halo_size)
                send = @views A[isend...]
                len = length(send)
                buf_view = @views sendbuf[(1:len) .+ nwrite]
                copyto!(buf_view, send)
                nwrite += len
            end
        end
        nwrites[i] = nwrite
    end

    KA.synchronize(backend)
end

function commhalo!(
    ::MultiProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    buf::MPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    comm = Topologies.comm()
    (; requests, sendbufs, recvbufs, nwrites) = buf

    neighbors = topo.neighbor_linear_ranks
    self = Topologies.linear_rank(topo)

    for (i, (neighbor, sendbuf, recvbuf, nwrite)) in enumerate(zip(neighbors, sendbufs, recvbufs, nwrites))
        if neighbor == self
            continue
        end

        send = @views sendbuf[1:nwrite]
        recv = @views recvbuf[1:nwrite]
        # Do not use reshape here, MultiRequest is not an array.
        lis = LinearIndices((N_COMMS, Topologies.N_DIRECTIONS))

        # Pair tags to coordinate send/recv directions
        # e.g., i=1 for x+, i=2 for x-; in general, i=2n-1 and i=2n form a reverse pair
        sendtag, recvtag = if i % 2 == 1
            i, -i
        else
            1 - i, i - 1
        end
        sendtag += OFFSET
        recvtag += OFFSET

        MPI.Isend(send, comm, requests[lis[1, i]]; dest=neighbor, tag=sendtag)
        MPI.Irecv!(recv, comm, requests[lis[2, i]]; source=neighbor, tag=recvtag)
    end
end

function waithalo(::MultiProcessing, buf::MPHaloBuffer)
    while !MPI.Testall(buf.requests)
        yield()
    end
end

function unpackhalo!(
    ::MultiProcessing,
    style::HaloStyle,
    As::Tuple{Vararg{AbstractArray{<:Real,N_DIMS}}},
    buf::MPHaloBuffer,
    dims::NTuple{N_DIMS,Integer},
    halo_size::Integer,
    backend::KA.Backend,
    topo::Topology,
)
    (; recvbufs, nwrites) = buf

    halo_size_comm = Topologies.halo_size_comm(style)
    directions = comm_directions(style)
    neighbors = topo.neighbor_linear_ranks
    self = Topologies.linear_rank(topo)

    for (i, (direction, neighbor, recvbuf)) in enumerate(zip(directions, neighbors, recvbufs))
        nread = 0
        for A in As
            if neighbor == self
                continue
            end
            irecv = offset_recv_indices(direction, dims, halo_size_comm, halo_size)
            recv = @views A[irecv...]
            len = length(recv)
            len + nread <= nwrites[i] ||
                throw(ArgumentError("Halo read size $(len + nread) exceeds write size $(nwrites[i])"))
            buf_view = @views recvbuf[(1:len) .+ nread]
            copyto!(recv, buf_view)
            nread += len
        end
    end

    KA.synchronize(backend)
end

function send_indices_to_other(direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, r::Integer)
    map(direction, dims) do e, n
        e == 0 && return 1:n
        e == +1 && return (n - r + 1):n
        e == -1 && return 1:r
        throw(ArgumentError("invalid direction $direction"))
    end
end

function send_indices_to_self(direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, r::Integer)
    map(direction, dims) do e, n
        e == 0 && return 1:n
        e == +1 && return 1:r
        e == -1 && return (n - r + 1):n
        throw(ArgumentError("invalid direction $direction"))
    end
end

function recv_indices(direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, r::Integer)
    map(direction, dims) do e, n
        e == 0 && return 1:n
        e == +1 && return (n + 1):(n + r)
        e == -1 && return (1 - r):0
        throw(ArgumentError("invalid direction $direction"))
    end
end

function offset_indices(indices::Tuple{Vararg{AbstractRange}}, r::Integer)
    map(indices) do range
        range .+ r
    end
end

function offset_send_indices_to_self(
    direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, halo_size_comm::Integer, halo_size_full::Integer
)
    offset_indices(send_indices_to_self(direction, dims, halo_size_comm), halo_size_full)
end

function offset_send_indices_to_other(
    direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, halo_size_comm::Integer, halo_size_full::Integer
)
    offset_indices(send_indices_to_other(direction, dims, halo_size_comm), halo_size_full)
end

function offset_recv_indices(
    direction::Tuple{Vararg{Integer}}, dims::Tuple{Vararg{Integer}}, halo_size_comm::Integer, halo_size_full::Integer
)
    offset_indices(recv_indices(direction, dims, halo_size_comm), halo_size_full)
end

function faces(domain::NTuple{N_DIMS,AbstractUnitRange}, core::NTuple{N_DIMS,AbstractUnitRange})
    range_x_domain, range_y_domain, range_z_domain = domain
    range_x_core, range_y_core, range_z_core = core

    start_x_domain, start_y_domain, start_z_domain = first.(domain)
    start_x_core, start_y_core, start_z_core = first.(core)
    stop_x_domain, stop_y_domain, stop_z_domain = last.(domain)
    stop_x_core, stop_y_core, stop_z_core = last.(core)

    (
        (range_x_domain, range_y_domain, start_z_domain:(start_z_core - 1)),
        (range_x_domain, range_y_domain, (stop_z_core + 1):stop_z_domain),
        (range_x_domain, start_y_domain:(start_y_core - 1), range_z_core),
        (range_x_domain, (stop_y_core + 1):stop_y_domain, range_z_core),
        (start_x_domain:(start_x_core - 1), range_y_core, range_z_core),
        ((stop_x_core + 1):stop_x_domain, range_y_core, range_z_core),
    )
end
