#=
  *--  --*--  --*--  --*
 /| 26  /| 14  /| 21  /|
*--  --*--  --*--  --*
/| 08  /| 05  /| 07  /|
*--  --*--  --*--  --*
/| 23  /| 11  /| 19  /|
*--  --*--  --*--  --*
|      |      |      |

  |      |      |      |
  *--  --*--  --*--  --*       z+ |  / x-
|/| 17 |/| 02 |/| 16 |/|          | /
*--  --*--  --*--  --*      y-    |/
|/| 04 |/| 00 |/| 03 |/|      ------*------
*--  --*--  --*--  --*             /|    y+
|/| 18 |/| 01 |/| 15 |/|            / |
*--  --*--  --*--  --*          x+ /  | z-
|      |      |      |

  |      |      |      |
  *--  --*--  --*--  --*
|/  20 |/  13 |/  24 |/
*--  --*--  --*--  --*
|/  09 |/  06 |/  10 |/
*--  --*--  --*--  --*
|/  22 |/  12 |/  25 |/
*--  --*--  --*--  --*
=#

"""
    DIRECTIONS

26 neighbor directions (6 face, 12 edge, 8 vertex) for 3D halo exchange.
"""
const DIRECTIONS = (
    # Face
    (+1, 0, 0),
    (-1, 0, 0),
    (0, +1, 0),
    (0, -1, 0),
    (0, 0, +1),
    (0, 0, -1),
    # Edge
    (0, +1, +1),
    (0, -1, -1),
    (0, +1, -1),
    (0, -1, +1),
    (+1, 0, +1),
    (-1, 0, -1),
    (-1, 0, +1),
    (+1, 0, -1),
    (+1, +1, 0),
    (-1, -1, 0),
    (+1, -1, 0),
    (-1, +1, 0),
    # Vertex
    (+1, +1, +1),
    (-1, -1, -1),
    (-1, +1, +1),
    (+1, -1, -1),
    (+1, -1, +1),
    (-1, +1, -1),
    (+1, +1, -1),
    (-1, -1, +1),
)
const N_DIRECTIONS = length(DIRECTIONS)
const N_DIMS = 3

"""
    Topology

3D Cartesian process topology with periodic boundaries.
"""
struct Topology{MPX,MPY,MPZ,NP}
    linear_rank           :: Int
    linear_ranks          :: NTuple{MPX,NTuple{MPY,NTuple{MPZ,Int}}}
    neighbor_linear_ranks :: NTuple{N_DIRECTIONS,Int}
    cart_rank             :: NTuple{N_DIMS,Int}
    cart_ranks            :: NTuple{NP,NTuple{N_DIMS,Int}}
end

"""
    Topology(proc_dims; reorder=true, multiprocessing=MPI.Initialized())

Map processes to 3D grid and provide neighbor ranks for halo exchange.
When `proc_dims = [0,0,0]`, MPI auto-determines dimensions. MPI must be
initialized before creating a topology in multiprocessing mode.
"""
function Topology(proc_dims::Vector{<:Integer}=[0, 0, 0]; reorder::Bool=true, multiprocessing::Bool=MPI.Initialized())
    if multiprocessing
        _topology_multiprocessing(proc_dims, reorder)
    else
        _topology_singleprocessing(proc_dims)
    end
end

"""
    comm()

MPI communicator for topology operations (requires MPI.Initialized()).
"""
function comm()
    MPI.Initialized() || throw(ArgumentError("MPI must be initialized"))
    MPI.COMM_WORLD
end

"""
    isroot(topo)

Check if current process is root (rank 0).
"""
function isroot(topo::Topology)
    linear_rank(topo) == 0
end

"""
    device!(backend; multiprocessing=MPI.Initialized())

Configure MPI-aware device selection for backend.
"""
function device!(backend::KA.Backend; multiprocessing::Bool=MPI.Initialized())
    if multiprocessing
        comm = Topologies.comm()
        linear_rank = MPI.Comm_rank(comm)

        select_device!(backend, comm, linear_rank)
    end
end

function select_device!(::CPU, comm::MPI.Comm, rank::Integer) end

function select_device!(backend::GPU, comm::MPI.Comm, rank::Integer)
    MPI.Initialized() || throw(ArgumentError("MPI must be initialized before selecting device, call MPI.Init() first"))

    lgrid = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank)
    lrank = MPI.Comm_rank(lgrid)
    if lrank < KA.ndevices(backend)
        KA.device!(backend, lrank + 1)
    end
end

"""
    proc_dims(topo)

Get process dimensions.
"""
function proc_dims(::Topology{MPX,MPY,MPZ}) where {MPX,MPY,MPZ}
    MPX, MPY, MPZ
end

"""
    proc_size(topo)

Get total number of processes.
"""
function proc_size(topo::Topology)
    prod(proc_dims(topo))
end

"""
    linear_rank(topo)

Get linear rank of current process.
"""
function linear_rank(topo::Topology)
    topo.linear_rank
end

"""
    linear_ranks(topo)

Iterable of all linear ranks (0 to proc_size-1).
"""
function linear_ranks(topo::Topology)
    0:(proc_size(topo) - 1)
end

"""
    cart_rank(topo)

Get Cartesian coordinates of current process.
"""
function cart_rank(topo::Topology)
    topo.cart_rank
end

"""
    cart_to_linear_rank(topo, cart_rank)

Convert Cartesian coordinates to linear rank.
"""
function cart_to_linear_rank(topo::Topology, cart_rank::NTuple{N_DIMS,Integer})
    ipx, ipy, ipz = cart_rank
    topo.linear_ranks[ipx + 1][ipy + 1][ipz + 1]
end

"""
    linear_to_cart_rank(topo, linear_rank)

Convert linear rank to Cartesian coordinates.
"""
function linear_to_cart_rank(topo::Topology, linear_rank::Integer)
    topo.cart_ranks[linear_rank + 1]
end

"""
    each_cart_rank(topo)

Get iterator over all Cartesian ranks.
"""
function each_cart_rank(topo::Topology)
    proc_dims = Topologies.proc_dims(topo)
    Iterators.product(map(mp -> 0:(mp - 1), proc_dims)...)
end

function _topology_multiprocessing(proc_dims::Vector{<:Integer}, reorder::Bool)
    MPI.Initialized() || throw(ArgumentError("MPI must be initialized before creating topology, call MPI.Init() first"))

    comm = Topologies.comm()
    proc_size = MPI.Comm_size(comm)

    length(proc_dims) == N_DIMS || throw(ArgumentError("proc_dims must have length $N_DIMS"))

    if all(proc_dims .== 0)
        proc_dims = MPI.Dims_create(proc_size, proc_dims)
    else
        proc_size == prod(proc_dims) ||
            throw(ArgumentError("proc_dims product $(prod(proc_dims)) does not match process count $proc_size"))
    end

    periodic = ones(Bool, N_DIMS)
    cart_comm = MPI.Cart_create(comm, proc_dims; periodic, reorder)

    linear_rank = MPI.Comm_rank(cart_comm)
    cart_rank = MPI.Cart_coords(cart_comm, linear_rank)

    linear_ranks = generate_linear_ranks(cart_comm, proc_dims)
    neighbor_linear_ranks = generate_neighbors(cart_comm, cart_rank, proc_dims, periodic)
    cart_ranks = generate_cart_ranks(cart_comm, proc_size)

    cart_rank = Int.(Tuple(cart_rank))
    Topology(linear_rank, linear_ranks, neighbor_linear_ranks, cart_rank, cart_ranks)
end

function _topology_singleprocessing(proc_dims::Vector{<:Integer})
    all(proc_dims .== 0) ||
        all(proc_dims .== 1) ||
        throw(ArgumentError("In single-processing mode, proc_dims must be all zeros or ones."))
    Topology(0, (((0,),),), ntuple(_ -> 0, N_DIRECTIONS), (0, 0, 0), ((0, 0, 0),))
end

function generate_neighbors(
    comm_cart::MPI.Comm,
    coords::AbstractVector{<:Integer},
    proc_dims::AbstractVector{<:Integer},
    periodic::AbstractVector{Bool},
)
    cart = collect(coords)
    map(DIRECTIONS) do disp
        candidate = cart .+ disp
        out_of_bounds = any(i -> !periodic[i] && (candidate[i] < 0 || candidate[i] >= proc_dims[i]), 1:N_DIMS)
        out_of_bounds ? -1 : MPI.Cart_rank(comm_cart, candidate)
    end
end

function generate_linear_ranks(comm_cart::MPI.Comm, proc_dims::AbstractVector{<:Integer})
    ntuple(proc_dims[1]) do ipx
        ntuple(proc_dims[2]) do ipy
            ntuple(proc_dims[3]) do ipz
                MPI.Cart_rank(comm_cart, [ipx - 1, ipy - 1, ipz - 1])
            end
        end
    end
end

function generate_cart_ranks(comm_cart::MPI.Comm, proc_length::Integer)
    ntuple(proc_length) do linear_rank
        Tuple(MPI.Cart_coords(comm_cart, linear_rank - 1))
    end
end
