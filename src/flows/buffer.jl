Utils.@concretize struct RCFBuffer{T<:Real}
    #! format: off
    Us_coarse_l :: AbstractArray{T,4}
    Us_coarse_g :: Tuple3{Array{T,3}}
    Us_hat      :: Tuple3{Array{Complex{T},3}}
    plan        :: FFTW.rFFTWPlan
    sendbuf     :: Union{Array{T,4},Nothing}
    recvbuf     :: Union{Vector{T},Nothing}
    temp        :: Array{Complex{T},3}
    #! format: on
end

function RCFBuffer(forcing::RCForcing, backend::KA.Backend, topo::Topologies.Topology)
    coarse_dims_g = Flows.coarse_dims_g(forcing)
    coarse_dims_l = Flows.coarse_dims_l(forcing, topo)

    Us_coarse_l = KA.zeros(backend, LCS.FP, coarse_dims_l..., LCS.N_DIMS)
    Us_coarse_g = Utils.@ntuple(zeros(LCS.FP, coarse_dims_g...), LCS.N_DIMS)
    plan = FFTW.plan_rfft(first(Us_coarse_g))
    Us_hat = Ref(plan) .* Us_coarse_g
    temp = FT.alloc_for_crop(Complex{LCS.FP}, forcing.kmax)

    if Topologies.is_multi_processing(Topologies.processing(topo))
        sendbuf = zeros(LCS.FP, coarse_dims_l..., LCS.N_DIMS)
        recvbuf = zeros(LCS.FP, prod(coarse_dims_g) * LCS.N_DIMS)
        RCFBuffer(Us_coarse_l, Us_coarse_g, Us_hat, plan, sendbuf, recvbuf, temp)
    else
        RCFBuffer(Us_coarse_l, Us_coarse_g, Us_hat, plan, nothing, nothing, temp)
    end
end

ForcingBuffer(::Forcing, config::LCS.Config, topo::Topologies.Topology) = nothing

function ForcingBuffer(forcing::RCForcing, config::LCS.Config, topo::Topologies.Topology)
    RCFBuffer(forcing, config.backend, topo)
end

Utils.@concretize struct Fields
    Us    :: Tuple3{Field}
    Us2   :: Tuple3{Field}
    P     :: Field
    dUdts :: Tuple3{Field}
end

function Fields(config::LCS.Config, topo::Topologies.Topology)
    Fields(config.backend, config.grid, topo)
end
function Fields(backend::KA.Backend, grid::LCS.Grid, topo::Topologies.Topology)
    dims = LCS.dims_l(grid.dims, topo)
    dims_with_halo = map(n -> n + 2 * grid.halo_size, dims)

    Us = Utils.@ntuple(KA.zeros(backend, LCS.FP, dims_with_halo...), LCS.N_DIMS)
    Us2 = Utils.@ntuple(KA.zeros(backend, LCS.FP, dims_with_halo...), LCS.N_DIMS)
    P = KA.zeros(backend, LCS.FP, dims_with_halo...)
    dUdts = Utils.@ntuple(KA.zeros(backend, LCS.FP, dims...), LCS.N_DIMS)

    Fields(Us, Us2, P, dUdts)
end

const Scratches{N} = NTuple{N,Field}
const N_SCRATCHES = 6

function Scratches(config::LCS.Config, topo::Topologies.Topology)
    Scratches(config.backend, config.grid, topo)
end
function Scratches(backend::KA.Backend, grid::LCS.Grid, topo::Topologies.Topology)
    dims = LCS.dims_l(grid.dims, topo)
    dims_with_halo = map(n -> n + 2 * grid.halo_size, dims)

    Utils.@ntuple(KA.zeros(backend, LCS.FP, dims_with_halo...), N_SCRATCHES)
end

Utils.@concretize struct IntegralLengthBuffer
    #! format: off
    R        :: AbstractArray{<:Real,3}
    Ls       :: AbstractArray{<:Real,2}
    U_hat    :: AbstractArray{<:Complex,3}
    plan     :: AbstractFFTs.Plan
    #! format: on
end

function IntegralLengthBuffer(config::LCS.Config, topo::Topologies.Topology)
    config.flow.stat.integral_length || return nothing
    IntegralLengthBuffer(config.backend, config.grid, topo)
end
function IntegralLengthBuffer(backend::KA.Backend, grid::LCS.Grid, topo::Topologies.Topology)
    n1, ns... = dims = LCS.dims_l(grid, topo)
    R = KA.zeros(backend, LCS.FP, dims...)
    Ls = similar(R, ns...)
    plan = FFTW.plan_rfft(R, 1)
    U_hat = plan * R
    IntegralLengthBuffer(R, Ls, U_hat, plan)
end

Utils.@concretize struct StatBuffer
    integral_length :: Union{IntegralLengthBuffer,Nothing}
end

function StatBuffer(config::LCS.Config, topo::Topologies.Topology)
    StatBuffer(IntegralLengthBuffer(config, topo))
end

"""
    FlowBuffer <: LCS.AbstractBuffer

Main buffer struct for flow computations.
"""
Utils.@concretize struct FlowBuffer <: LCS.AbstractBuffer
    fields    :: Fields
    scratches :: Scratches
    forcing   :: Union{RCFBuffer,Nothing}
    stat      :: StatBuffer
end

function FlowBuffer(config::LCS.Config, topo::Topologies.Topology)
    fields = Fields(config, topo)
    scratches = Scratches(config, topo)
    forcingbuf = ForcingBuffer(config.flow.forcing, config, topo)
    statbuf = StatBuffer(config, topo)
    FlowBuffer(fields, scratches, forcingbuf, statbuf)
end
