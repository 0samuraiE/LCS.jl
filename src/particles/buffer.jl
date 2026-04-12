Utils.@concretize struct Properties
    ids    :: Property{<:Integer}
    diams  :: Property{<:Real}
    xss    :: Tuple3{<:Property{<:Real}}
    xss2   :: Tuple3{<:Property{<:Real}}
    uss    :: Tuple3{<:Property{<:Real}}
    uss2   :: Tuple3{<:Property{<:Real}}
    dudtss :: Tuple3{<:Property{<:Real}}
end

Adapt.@adapt_structure Properties

function Properties(config::LCS.Config, topo::Topologies.Topology, iprofile::Integer=1)
    Properties(config.backend, config.particles[iprofile].population, topo)
end
function Properties(backend::KA.Backend, population::Population, topo::Topologies.Topology)
    ntotal = ntotal_l(population, topo)

    ids = KA.zeros(backend, LCS.IP, ntotal)
    fill!(ids, INVALID)

    diams = KA.zeros(backend, LCS.FP, ntotal)
    xss = LCS.@ntuple(KA.zeros(backend, LCS.FP, ntotal), LCS.N_DIMS)
    xss2 = LCS.@ntuple(KA.zeros(backend, LCS.FP, ntotal), LCS.N_DIMS)
    uss = LCS.@ntuple(KA.zeros(backend, LCS.FP, ntotal), LCS.N_DIMS)
    uss2 = LCS.@ntuple(KA.zeros(backend, LCS.FP, ntotal), LCS.N_DIMS)
    dudtss = LCS.@ntuple(KA.zeros(backend, LCS.FP, ntotal), LCS.N_DIMS)

    Properties(ids, diams, xss, xss2, uss, uss2, dudtss)
end

function Base.Tuple(props::Properties)
    (; ids, diams, xss, xss2, uss, uss2, dudtss) = props
    (ids, diams, xss..., xss2..., uss..., uss2..., dudtss...)
end

const N_PROPERTIES = 17

function Utils.sortbyperm!(props::Properties, perm::Property{<:Integer}, copy::Property{<:Real}, nvalid::Integer)
    for A in Tuple(props)
        Utils.sortbyperm!(A, perm, copy, nvalid)
    end
end

Utils.@concretize struct CellIndexBuffer
    perm   :: Property{<:Integer}
    hashes :: Property{<:Integer}
    starts :: AbstractVector{<:Integer}
    stops  :: AbstractVector{<:Integer}
    copy   :: Property{<:Real}
    hasher :: Utils.Hasher{LCS.N_DIMS}
end

Adapt.@adapt_structure CellIndexBuffer

function CellIndexBuffer(config::LCS.Config, topo::Topologies.Topology, iprofile::Integer=1)
    CellIndexBuffer(
        config.backend, config.grid, config.particles[iprofile].population, config.particles[iprofile].cell, topo
    )
end
function CellIndexBuffer(
    backend::KA.Backend, grid::LCS.Grid, population::Population, cell::Cell, topo::Topologies.Topology
)
    n_total = ntotal_l(population, topo)
    cell_dims = LCS.dims_l(cell, grid, topo)

    perm = KA.zeros(backend, LCS.IP, n_total)
    hashes = KA.zeros(backend, LCS.IP, n_total)
    starts = KA.zeros(backend, LCS.IP, prod(cell_dims .+ 2))
    stops = KA.zeros(backend, LCS.IP, prod(cell_dims .+ 2))
    copy = KA.zeros(backend, LCS.FP, n_total)
    hasher = Utils.Hasher(map(n -> 0:(n + 1), cell_dims))

    CellIndexBuffer(perm, hashes, starts, stops, copy, hasher)
end

const N_ONE_SIDES = 2
const N_ALL_SIDES = 6

Utils.@concretize struct CommBuffer
    indices :: Property{<:Integer}
    mask    :: Property{Bool}
    scan    :: Property{<:Integer}
    perm    :: Property{<:Integer}
    copy    :: Property{<:Real}
    sends   :: NTuple{N_ONE_SIDES,<:AbstractVector{<:Real}}
    recvs   :: NTuple{N_ONE_SIDES,<:AbstractVector{<:Real}}
end

Adapt.@adapt_structure CommBuffer

function CommBuffer(config::LCS.Config, topo::Topologies.Topology, cibuf::CellIndexBuffer, iprofile::Integer=1)
    CommBuffer(config.backend, config.particles[iprofile].population, topo, cibuf)
end
function CommBuffer(
    backend::KA.Backend,
    population::Population,
    topo::Topologies.Topology,
    cibuf::Union{CellIndexBuffer,Nothing}=nothing,
)
    if Topologies.is_multi_processing(topo)
        ntotal = ntotal_l(population, topo)
        nexchange = nexchange_l(population, topo)

        indices = isnothing(cibuf) ? KA.zeros(backend, LCS.IP, ntotal) : cibuf.hashes
        mask = KA.zeros(backend, Bool, ntotal)
        scan = KA.zeros(backend, LCS.IP, ntotal)
        perm = isnothing(cibuf) ? KA.zeros(backend, LCS.IP, ntotal) : cibuf.perm
        copy = isnothing(cibuf) ? KA.zeros(backend, LCS.FP, ntotal) : cibuf.copy
        sends = Utils.@ntuple(KA.zeros(backend, LCS.FP, nexchange * N_PROPERTIES), N_ONE_SIDES)
        recvs = Utils.@ntuple(KA.zeros(backend, LCS.FP, nexchange * N_PROPERTIES), N_ONE_SIDES)

        CommBuffer(indices, mask, scan, perm, copy, sends, recvs)
    else
        nothing
    end
end

"""
    ParticleBuffer <: LCS.AbstractBuffer

Main buffer struct for particle computations.
"""
Utils.@concretize struct ParticleBuffer <: LCS.AbstractBuffer
    props :: Properties
    ci    :: CellIndexBuffer
    comm  :: Union{CommBuffer,Nothing}
end

function ParticleBuffer(config::LCS.Config, topo::Topologies.Topology, iprofile::Integer=1)
    props = Properties(config, topo, iprofile)
    cibuf = CellIndexBuffer(config, topo, iprofile)
    commbuf = CommBuffer(config, topo, cibuf, iprofile)
    ParticleBuffer(props, cibuf, commbuf)
end
