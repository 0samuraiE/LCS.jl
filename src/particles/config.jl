"""
    Population

Particle population configuration.

# Fields
- `valid`: Total number of valid particles across all processes
- `invalid_ratio`: Ratio of buffer space for invalid particles
- `exchange_ratio`: Ratio of buffer space for particle exchange between processes
"""
struct Population
    valid          :: Int
    invalid_ratio  :: Float64
    exchange_ratio :: Float64
end
@composite Population

nvalid_g(p::Population) = p.valid

nvalid_l(p::Population, topo::Topologies.Topology) = div(nvalid_g(p), Topologies.proc_size(topo))

ntotal_l(p::Population, topo::Topologies.Topology) = ceil(Int, nvalid_l(p, topo) * (1 + p.invalid_ratio))

nexchange_l(p::Population, topo::Topologies.Topology) = ceil(Int, nvalid_l(p, topo) * (1 + p.exchange_ratio))

"""
    ParticleParams

Physical parameters for particle dynamics.

# Fields
- `gravity`: Gravitational acceleration (non-dimensional)
- `density_ratio`: Particle to fluid density ratio (ρₚ/ρf)
"""
struct ParticleParams
    gravity       :: Float64
    density_ratio :: Float64
end
@composite ParticleParams

taup(pp::ParticleParams, fp::Flows.FlowParams, diam::Real) = 1 / 18 * pp.density_ratio * Flows.Re(fp) * diam^2

"""
    Cell

Cell-based particle tracking configuration.

# Fields
- `grid_ratio`: Coarsening ratio relative to flow grid
"""
struct Cell
    grid_ratio :: Int
end
@composite Cell

LCS.dims_g(c::Cell, grid::LCS.Grid) = div.(LCS.dims_g(grid), c.grid_ratio)

LCS.dims_l(c::Cell, grid::LCS.Grid, topo::Topologies.Topology) = div.(LCS.dims_l(grid, topo), c.grid_ratio)

LCS.spacings(c::Cell, grid::LCS.Grid) = LCS.spacings(grid) .* c.grid_ratio

function cell_in_domain(ci::Integer3, cell::Cell, grid::LCS.Grid, topo::Topologies.Topology)
    all(1 .<= ci .<= LCS.dims_l(cell, grid, topo))
end

"""
    RDF

Radial distribution function configuration.

# Fields
- `contact_tolerance`: Tolerance for contact distance detection
- `bin_size`: Number of grid cells per bin
- `min_radius`: Minimum radius for binning
- `search_cell_size`: Number of cells to search around each particle
"""
struct RDF
    contact_tolerance :: Float64
    bin_size          :: Int
    min_radius        :: Float64
    search_cell_size  :: Int
end
@composite RDF

"""
    Spectrum

Particle number-density spectrum configuration.
"""
struct Spectrum
    max_wavenumber :: Int
    nsamples       :: Int
end
@composite Spectrum

"""
    Stat

Particle statistics configuration.

# Fields
- `rdf`: Radial distribution function configuration
- `spectrum`: Particle number-density spectrum configuration
"""
struct Stat
    rdf      :: RDF
    spectrum :: Spectrum
end
@composite Stat

"""
    DragModel

Abstract type for drag force models.
"""
abstract type DragModel end
@variant DragModel

"""
    LinearDrag <: DragModel

Linear Stokes drag model.
"""
struct LinearDrag <: DragModel end
@composite LinearDrag
@kind DragModel "linear" LinearDrag

coeff(::LinearDrag, Rep::Real) = 1.0

"""
    NonlinearDrag <: DragModel

Nonlinear drag model with Reynolds number correction using Schiller-Naumann correlation.
"""
struct NonlinearDrag <: DragModel end
@composite NonlinearDrag
@kind DragModel "nonlinear" NonlinearDrag

coeff(::NonlinearDrag, Rep::Real) = 1.0 + 0.15 * Rep^0.687

"""
    ParticleRestart

Particle restart configuration specifying source file and profile index.

# Fields
- `file`: Path to restart HDF5 file
- `iprofile`: Particle profile index (1-based)
"""
struct ParticleRestart
    file     :: String
    iprofile :: Int
end
@composite ParticleRestart

"""
    InitId

Abstract type for particle ID initialization methods.
"""
abstract type InitId end
@variant Union{InitId,ParticleRestart}
@kind Union{InitId,ParticleRestart} "restart" ParticleRestart

"""
    GenerateId <: InitId

Generate unique particle IDs for each particle.
"""
struct GenerateId <: InitId end
@composite GenerateId
@kind Union{InitId,ParticleRestart} "generate" GenerateId

"""
    InitPosition

Abstract type for particle position initialization methods.
"""
abstract type InitPosition end
@variant Union{InitPosition,ParticleRestart}
@kind Union{InitPosition,ParticleRestart} "restart" ParticleRestart

"""
    RandomPosition <: InitPosition

Initialize particle positions randomly in domain.
"""
struct RandomPosition <: InitPosition end
@composite RandomPosition
@kind Union{InitPosition,ParticleRestart} "random" RandomPosition

"""
    InitVelocity

Abstract type for particle velocity initialization methods.
"""
abstract type InitVelocity end
@variant Union{InitVelocity,ParticleRestart}
@kind Union{InitVelocity,ParticleRestart} "restart" ParticleRestart

"""
    RestVelocity <: InitVelocity

Initialize particles at rest (zero velocity).
"""
struct RestVelocity <: InitVelocity end
@composite RestVelocity
@kind Union{InitVelocity,ParticleRestart} "rest" RestVelocity

"""
    InitSize

Abstract type for particle size initialization methods.
"""
abstract type InitSize end
@variant Union{InitSize,ParticleRestart}
@kind Union{InitSize,ParticleRestart} "restart" ParticleRestart

"""
    ConstSize <: InitSize

Initialize all particles with constant diameter.

# Fields
- `diam_m`: Particle diameter [m]
"""
struct ConstSize <: InitSize
    diam_m :: Float64
end
@composite ConstSize
@kind Union{InitSize,ParticleRestart} "const" ConstSize

diam(cs::ConstSize, fparams::Flows.FlowParams) = cs.diam_m / fparams.L0_m

"""
    Init

Particle initialization configuration combining ID, position, velocity, and size initialization.

# Fields
- `id`: ID initialization method
- `position`: Position initialization method
- `velocity`: Velocity initialization method
- `size`: Size initialization method
"""
struct Init{
    InitIdT<:Union{InitId,ParticleRestart},
    InitPositionT<:Union{InitPosition,ParticleRestart},
    InitVelocityT<:Union{InitVelocity,ParticleRestart},
    InitSizeT<:Union{InitSize,ParticleRestart},
}
    id       :: InitIdT
    position :: InitPositionT
    velocity :: InitVelocityT
    size     :: InitSizeT
end
@composite Init

"""
    ParticleConfig <: LCS.AbstractConfig

Main configuration struct for particle simulation.

# Fields
- `population`: Particle population configuration
- `params`: Physical parameters
- `cell`: Cell-based indexing configuration
- `stat`: Statistics configuration
- `drag`: Drag force model
- `init`: Initialization configuration
"""
struct ParticleConfig{
    PopulationT<:Population,ParamsT<:ParticleParams,CellT<:Cell,StatT<:Stat,DragModelT<:DragModel,InitT<:Init
} <: LCS.AbstractConfig
    population :: PopulationT
    params     :: ParamsT
    cell       :: CellT
    stat       :: StatT
    drag       :: DragModelT
    init       :: InitT
end
@composite ParticleConfig
