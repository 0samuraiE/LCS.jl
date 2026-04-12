"""
    FlowParams

Flow parameters defining the physical properties of the simulation.

# Fields
- `U0_m_s`: Reference velocity scale [m/s]
- `L0_m`: Reference length scale [m]
- `NU_m2_s`: Kinematic viscosity [m²/s]
"""
struct FlowParams
    U0_m_s  :: Float64
    L0_m    :: Float64
    NU_m2_s :: Float64
end
@composite FlowParams

Re(params::FlowParams) = params.U0_m_s * params.L0_m / params.NU_m2_s

"""
    Init

Abstract type for flow field initialization methods.
"""
abstract type Init end
@variant Union{Init,LCS.Restart}
@kind Union{Init,LCS.Restart} "restart" LCS.Restart

"""
    RandomFlow <: Init

Initialize flow field with random velocity fluctuations.

# Fields
- `seed`: Random seed for reproducible initialization
"""
struct RandomFlow <: Init
    seed :: Int
end
@composite RandomFlow
@kind Union{Init,LCS.Restart} "random-flow" RandomFlow

"""
    IdealFlow <: Init

Initialize flow field with ideal (analytical) velocity solution.
"""
struct IdealFlow <: Init end
@composite IdealFlow
@kind Union{Init,LCS.Restart} "ideal-flow" IdealFlow

"""
    Upsample <: Init

Initialize flow field by upsampling from coarser resolution restart file.

# Fields
- `file`: Path to restart file to upsample from
"""
struct Upsample <: Init
    file :: String
end
@composite Upsample
@kind Union{Init,LCS.Restart} "upsample" Upsample

"""
    Couple

Abstract type for pressure-velocity coupling methods.
"""
abstract type Couple end
@variant Couple

"""
    HSMAC <: Couple

HSMAC (Highly Simplified Marker and Cell) pressure-velocity coupling method.

# Fields
- `time_blocking`: Enable time-blocking
- `overlap`: Enable computation-communication overlap
- `coloring`: Coloring scheme for parallel iteration
- `omega`: SOR relaxation parameter (1.0 < ω < 2.0)
- `epsp`: Pressure convergence tolerance
- `itrmin`: Minimum iteration count
- `itrmax`: Maximum iteration count
"""
struct HSMAC{ColoringT<:Parallel.Coloring} <: Couple
    time_blocking :: Bool
    overlap       :: Bool
    coloring      :: ColoringT
    omega         :: Float64
    epsp          :: Float64
    itrmin        :: Int
    itrmax        :: Int
end
@composite HSMAC
@kind Couple "hsmac" HSMAC

"""
    Forcing

Abstract type for flow forcing methods.
"""
abstract type Forcing end
@variant Forcing

"""
    NoForcing <: Forcing

No external forcing applied to the flow.
"""
struct NoForcing <: Forcing end
@composite NoForcing
@kind Forcing "no" NoForcing

"""
    LinearForcing

Abstract type for linear forcing methods.
"""
abstract type LinearForcing <: Forcing end

"""
    ConstantPowerLF <: LinearForcing

Linear forcing that maintains constant power injection into the flow.

# Fields
- `power`: Target power injection rate
"""
struct ConstantPowerLF <: LinearForcing
    power :: Float64
end
@composite ConstantPowerLF
@kind Forcing "linear-power" ConstantPowerLF

"""
    RCForcing

Reduced Communication Forcing method for maintaining turbulent flow.
"""
abstract type RCForcing <: Forcing end

coarse_dims_g(forcing::RCForcing) = Utils.@ntuple(forcing.coarse_size, LCS.N_DIMS)
coarse_dims_l(forcing::RCForcing, topo::Topologies.Topology) = LCS.dims_l(coarse_dims_g(forcing), topo)
filter_size(forcing::RCForcing, grid::LCS.Grid) = div(first(LCS.dims_g(grid)), forcing.coarse_size)
filter_dims(forcing::RCForcing, grid::LCS.Grid) = div.(LCS.dims_g(grid), coarse_dims_g(forcing))
phase_shift(forcing::RCForcing, grid::LCS.Grid) = filter_size(forcing, grid) / 2

"""
    EnergyPreserveRCF <: RCForcing

Reduced Communication Forcing that maintains constant kinetic energy.

# Fields
- `energy`: Time-varying target energy
- `kmin`: Minimum wavenumber shell for forcing
- `kmax`: Maximum wavenumber shell for forcing
- `coarse_size`: Size of coarse grid for forcing generation
"""
struct EnergyPreserveRCF <: RCForcing
    energy      :: Float64
    kmin        :: Int
    kmax        :: Int
    coarse_size :: Int

    function EnergyPreserveRCF(energy::Float64, kmin::Int, kmax::Int, coarse_size::Int)
        kmin >= 1 || throw(ArgumentError("kmin must be at least 1, got $kmin"))
        kmax >= kmin || throw(ArgumentError("kmax must be at least kmin, got kmax=$kmax kmin=$kmin"))
        new(energy, kmin, kmax, coarse_size)
    end
end
@composite EnergyPreserveRCF
@kind Forcing "rc-energy" EnergyPreserveRCF

"""
    FlowStatConfig

Configuration for flow statistics computation.
"""
struct Stat
    integral_length :: Bool
end
@composite Stat

"""
    FlowConfig <: LCS.AbstractConfig

Main configuration struct for flow simulation.

# Fields
- `params`: Flow physical parameters
- `init`: Flow field initialization method
- `couple`: Pressure-velocity coupling method
- `forcing`: Flow forcing method
- `stat`: Flow statistics options
"""
struct FlowConfig <: LCS.AbstractConfig
    params  :: FlowParams
    init    :: Union{Init,LCS.Restart}
    couple  :: Couple
    forcing :: Forcing
    stat    :: Stat
end
@composite FlowConfig
