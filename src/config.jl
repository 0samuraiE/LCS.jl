"""
    Mode

Abstract base type for simulation modes.
"""
abstract type Mode end
@variant Mode

"""
    FlowMode <: Mode

Flow-only simulation mode (no particles).
"""
struct FlowMode <: Mode end
@composite FlowMode
@kind Mode "flow" FlowMode

"""
    FlowParticleMode <: Mode

Flow-particle coupled simulation mode.
"""
struct FlowParticleMode <: Mode end
@composite FlowParticleMode
@kind Mode "flow-particle" FlowParticleMode

"""
    TopoDraft

Process topology configuration.

# Fields
- `proc_dims`: Process grid dimensions
- `reorder`: Enable process reordering
"""
struct TopoDraft
    proc_dims :: Tuple3{Int}
    reorder   :: Bool
end
@composite TopoDraft
function PolySerde.deserialize(::PolySerde.Composite, ::Type{TopoDraft}, x::AbstractDict{String})
    passthrough = Dict{String,Any}("proc_dims" => Tuple3{Int}(x["proc_dims"]))
    PolySerde.deserialize(PolySerde.Composite(), TopoDraft, x, passthrough)
end
function PolySerde.serialize(::PolySerde.Composite, ::Type{TopoDraft}, x::TopoDraft, dicttype::PolySerde.SerdeDictType)
    passthrough = Dict{String,Any}("proc_dims" => collect(x.proc_dims))
    PolySerde.serialize(PolySerde.Composite(), TopoDraft, x, dicttype, passthrough)
end

"""
    Grid

Spatial grid configuration.

# Fields
- `dims`: Grid dimensions
- `halo_size`: Halo region size
"""
struct Grid
    dims      :: Tuple3{Int}
    halo_size :: Int
end
@composite Grid
function PolySerde.deserialize(::PolySerde.Composite, ::Type{Grid}, x::AbstractDict{String})
    passthrough = Dict{String,Any}("dims" => Tuple3{Int}(x["dims"]))
    PolySerde.deserialize(PolySerde.Composite(), Grid, x, passthrough)
end
function PolySerde.serialize(::PolySerde.Composite, ::Type{Grid}, x::Grid, dicttype::PolySerde.SerdeDictType)
    passthrough = Dict{String,Any}("dims" => collect(x.dims))
    PolySerde.serialize(PolySerde.Composite(), Grid, x, dicttype, passthrough)
end

Offsets.offset(g::Grid, i) = offset(g.halo_size, i)
Base.Broadcast.broadcastable(x::Grid) = Ref(x)

index_g(index_l::Integer3, dims_l::Integer3, cart_rank::Integer3) = @. index_l + cart_rank .* dims_l
function index_g(index_l::Integer3, dims_l::Integer3, topo::Topologies.Topology)
    index_g(index_l, dims_l, Topologies.cart_rank(topo))
end
# Do not implement this signature to avoid performance regression due to expensive divisions.
# index_g(index_l::Integer3, grid::Grid, topo::Topologies.Topology)

dims_g(grid::Grid) = grid.dims

dims_l(dims::Integer3, proc_dims::Integer3) = div.(dims, proc_dims)
dims_l(dims::Integer3, topo::Topologies.Topology) = dims_l(dims, Topologies.proc_dims(topo))
dims_l(grid::Grid, topo::Topologies.Topology) = dims_l(dims_g(grid), Topologies.proc_dims(topo))

spacings(grid::Grid) = LCS.DOMAIN_LENGTH ./ dims_g(grid)

lengths_g() = @ntuple(LCS.DOMAIN_LENGTH, LCS.N_DIMS)

lengths_l(lengths_g::Real3, proc_dims::Integer3) = lengths_g ./ proc_dims
lengths_l(topo::Topologies.Topology) = lengths_l(lengths_g(), Topologies.proc_dims(topo))

origins_g() = @ntuple(LCS.DOMAIN_ORIGIN, LCS.N_DIMS)

origins_l(origin_g::Real3, lengths_l::Real3, cart_rank::Integer3) = origin_g .+ lengths_l .* cart_rank
origins_l(topo::Topologies.Topology) = origins_l(origins_g(), lengths_l(topo), Topologies.cart_rank(topo))

ends_l(origin_l::Real3, lengths_l::Real3) = origin_l .+ lengths_l
ends_l(topo::Topologies.Topology) = ends_l(origins_l(topo), lengths_l(topo))

volume_g() = prod(lengths_g())

volume_l(topo::Topologies.Topology) = prod(lengths_l(topo))

indices_l(dim_l::Integer, rank::Integer) = (1:dim_l) .+ dim_l * rank
indices_l(dims_l::Integer3, cart_rank::Integer3) = indices_l.(dims_l, cart_rank)
indices_l(dim_l::Integer, topo::Topologies.Topology) = indices_l(dim_l, Topologies.linear_rank(topo))
indices_l(dims_l::Integer3, topo::Topologies.Topology) = indices_l(dims_l, Topologies.cart_rank(topo))
indices_l(grid::Grid, topo::Topologies.Topology) = indices_l(dims_l(grid, topo), Topologies.cart_rank(topo))

unhalo(A::AbstractArray, grid) = unhalo(A, grid.halo_size)
function unhalo(A::AbstractArray, halo_size::Integer)
    indices = map(n -> (halo_size + 1):(n - halo_size), size(A))
    view(A, indices...)
end

"""
    TimeStep

Time step parameters.

# Fields
- `cfl`: CFL number
- `dtmax`: Maximum time step
"""
struct TimeStep
    cfl   :: Float64
    dtmax :: Float64
end
@composite TimeStep

"""
    SimulateConfig

Simulation control parameters.

# Fields
- `last_step`: Final step number
- `interval_stat`: Statistics output interval
- `interval_restart`: Restart output interval
- `interval_gc`: GC interval
- `num_restart_files`: Restart files to retain
"""
struct SimulateConfig
    last_step         :: Int
    interval_stat     :: Int
    interval_restart  :: Int
    interval_gc       :: Int
    num_restart_files :: Int
end
@composite SimulateConfig

abstract type AbstractConfig end

"""
    Config

Main simulation configuration.

# Fields
- `mode`: Simulation mode
- `backend`: Compute backend
- `topo`: Process topology
- `grid`: Spatial grid
- `timestep`: Time step parameters
- `simulate`: Control parameters
- `flow`: Flow configuration
- `particle`: Particle configuration
"""
struct Config{
    ModeT<:Mode,
    BackendT<:KA.Backend,
    TopoDraftT<:TopoDraft,
    GridT<:Grid,
    TimeStepT<:TimeStep,
    SimulateT<:SimulateConfig,
    FlowT<:Union{AbstractConfig,Nothing},
    ParticlesT<:Union{Tuple{Vararg{AbstractConfig}},Nothing},
} <: AbstractConfig
    outdir    :: String
    resume    :: Bool
    mode      :: ModeT
    backend   :: BackendT
    topo      :: TopoDraftT
    grid      :: GridT
    timestep  :: TimeStepT
    simulate  :: SimulateT

    flow      :: FlowT
    particles :: ParticlesT
end
@composite Config

function PolySerde.deserialize(::PolySerde.Composite, ::Type{Config}, x::AbstractDict{String,Any})
    flow = if haskey(x, "flow") && !isnothing(x["flow"])
        PolySerde.deserialize(Flows.FlowConfig, x["flow"])
    else
        nothing
    end

    particles = if haskey(x, "particles") && !isnothing(x["particles"])
        ntuple(length(x["particles"])) do n
            PolySerde.deserialize(Particles.ParticleConfig, x["particles"][n])
        end
    else
        nothing
    end

    PolySerde.deserialize(PolySerde.Composite(), Config, x, Dict{String,Any}("flow" => flow, "particles" => particles))
end

function PolySerde.serialize(::PolySerde.Composite, ::Type{Config}, x::Config, dicttype::PolySerde.SerdeDictType)
    flow = if !isnothing(x.flow)
        PolySerde.serialize(Flows.FlowConfig, x.flow, dicttype)
    else
        nothing
    end

    particles = if !isnothing(x.particles)
        collect(
            ntuple(length(x.particles)) do n
                PolySerde.serialize(Particles.ParticleConfig, x.particles[n], dicttype)
            end,
        )
    else
        nothing
    end

    PolySerde.serialize(PolySerde.Composite(), Config, x, dicttype, dicttype("flow" => flow, "particles" => particles))
end

struct Restart
    file :: String
end
@composite Restart
