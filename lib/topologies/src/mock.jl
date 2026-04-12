module Mock
using ..Topologies

@inline cart_to_linear((x, y, z), (MPX, MPY, MPZ)) = x * MPY * MPZ + y * MPZ + z

@inline function linear_to_cart(r, (MPX, MPY, MPZ))
    z = mod(r, MPZ)
    y = mod(div(r, MPZ), MPY)
    x = div(r, MPY * MPZ)
    (x, y, z)
end

"""
    Topology(; linear_rank, proc_dims)

Create a mock topology for testing.
"""
function mocktopology(; linear_rank::Integer, proc_dims::NTuple{Topologies.N_DIMS,Integer})
    MPX, MPY, MPZ = proc_dims
    NP = MPX * MPY * MPZ

    lin2cart = lr -> linear_to_cart(lr, proc_dims)
    cart2lin = cr -> cart_to_linear(cr, proc_dims)

    cartesian_rank = lin2cart(linear_rank)

    linear_ranks = ntuple(ipx -> ntuple(ipy -> ntuple(ipz -> cart2lin((ipx - 1, ipy - 1, ipz - 1)), MPZ), MPY), MPX)

    neighbor_linear_ranks = ntuple(i -> begin
        dx, dy, dz = Topologies.DIRECTIONS[i]
        x = mod(cartesian_rank[1] + dx, MPX)
        y = mod(cartesian_rank[2] + dy, MPY)
        z = mod(cartesian_rank[3] + dz, MPZ)
        cart2lin((x, y, z))
    end, Topologies.N_DIRECTIONS)

    cartesian_ranks = ntuple(r -> lin2cart(r - 1), NP)

    Topologies.Topology(linear_rank, linear_ranks, neighbor_linear_ranks, cartesian_rank, cartesian_ranks)
end
end
