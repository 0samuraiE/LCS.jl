"""
    RDFStat

Radial distribution function statistics.

# Fields
- `gr_contact`: g(r) at contact distance
- `npairs_contact`: Number of pairs at contact
- `edges`: Bin edges
- `gr`: g(r) values per bin
- `npairs`: Number of pairs per bin
"""
@kwdef struct RDFStat
    #! format: off
    gr_contact     :: Float64         = 0.0
    npairs_contact :: Int             = 0
    edges          :: Vector{Float64} = Float64[]
    gr             :: Vector{Float64} = Float64[]
    npairs         :: Vector{Int}     = Int[]
    #! format: on
end
@composite RDFStat

PPrint.PrintStyle(::Type{<:RDFStat}) = PPrint.Tree()
PPrint.rename(::Val{:gr_contact}) = "g(r=R)"
PPrint.rename(::Val{:npairs_contact}) = "npair(r=R)"

"""
    ParticleStat <: LCS.AbstractStat

Particle statistics containing computed particle properties.

# Fields
- `x`: Position statistics in x-direction
- `y`: Position statistics in y-direction
- `z`: Position statistics in z-direction
- `u`: Velocity statistics in x-direction
- `v`: Velocity statistics in y-direction
- `w`: Velocity statistics in z-direction
- `diam`: Particle diameter statistics
- `τp`: Mean particle response time
- `St`: Stokes number (τp / τη)
- `rdf`: Radial distribution function statistics
"""
@kwdef struct ParticleStat <: LCS.AbstractStat
    #! format: off
    x    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    y    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    z    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    u    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    v    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    w    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    diam :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    τp   :: Float64              = 0.0
    St   :: Float64              = 0.0
    rdf  :: RDFStat              = RDFStat()
    #! format: on
end
@composite ParticleStat
PPrint.PrintStyle(::Type{<:ParticleStat}) = PPrint.Tree()

"""
    stat(fstat, buf, state, config, topo, iprofile)

Compute particle statistics.
"""
function stat(
    fstat::Flows.FlowStat,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    stat(;
        fstat,
        buf.particles[iprofile].props.xss,
        buf.particles[iprofile].props.uss,
        buf.particles[iprofile].props,
        cibuf=buf.particles[iprofile].ci,
        state.particles[iprofile].nvalid,
        config.grid,
        fparams=config.flow.params,
        pparams=config.particles[iprofile].params,
        config.particles[iprofile].population,
        config.particles[iprofile].cell,
        config.particles[iprofile].stat,
        config.backend,
        topo,
    )
end

function stat(;
    fstat::Flows.FlowStat,
    xss::Tuple3{Property{<:Real}},
    uss::Tuple3{Property{<:Real}},
    props::Properties,
    cibuf::CellIndexBuffer,
    nvalid::Integer,
    grid::LCS.Grid,
    fparams::Flows.FlowParams,
    pparams::ParticleParams,
    population::Population,
    cell::Cell,
    stat::Stat,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    xs, ys, zs = xss
    us, vs, ws = uss
    diams = props.diams

    @views begin
        x = LCS.Summary(xs[1:nvalid], topo)
        y = LCS.Summary(ys[1:nvalid], topo)
        z = LCS.Summary(zs[1:nvalid], topo)
        u = LCS.Summary(us[1:nvalid], topo)
        v = LCS.Summary(vs[1:nvalid], topo)
        w = LCS.Summary(ws[1:nvalid], topo)
        diam = LCS.Summary(diams[1:nvalid], topo)
    end

    τp = Particles.taup(pparams, fparams, diam.mean)
    St = τp / fstat.τη
    makeindex!(; xss, props, cibuf, nvalid, grid, cell, backend, topo)
    rdf = Particles.rdf_assume_ci(;
        diam=diam.mean, xss, props, cibuf, nvalid, grid, population, cell, topo, stat, backend
    )
    ParticleStat(; x, y, z, u, v, w, diam, τp, St, rdf)
end

"""
    density(buf, state, config, topo, iprofile)

Compute SPH (Smoothed Particle Hydrodynamics) density estimate for each particle.
"""
function density(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer)
    density(;
        buf.particles[iprofile].props.xss,
        buf.particles[iprofile].props,
        cibuf=buf.particles[iprofile].ci,
        state.particles[iprofile].nvalid,
        config.grid,
        config.particles[iprofile].cell,
        config.backend,
        topo,
    )
end

function density(;
    xss::Tuple3{Property{<:Real}},
    props::Properties,
    cibuf::CellIndexBuffer,
    nvalid::Integer,
    grid::LCS.Grid,
    cell::Cell,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    xs, ys, zs = xss
    (; hashes, starts, stops, hasher, copy) = cibuf

    densities = copy
    h = minimum(LCS.spacings(grid))

    @kernel function kernel!()
        tid = @index(Global, Linear)

        x1 = xs[tid], ys[tid], zs[tid]
        hash1 = hashes[tid]
        ci1 = Utils.decode(hasher, hash1)

        rho = 0.0

        if cell_in_domain(ci1, cell, grid, topo)
            for di in -1:1, dj in -1:1, dk in -1:1
                ci2 = ci1 .+ (di, dj, dk)
                hash2 = Utils.encode(hasher, ci2)

                for i2 in starts[hash2]:stops[hash2]
                    x2 = xs[i2], ys[i2], zs[i2]
                    d = Utils.norm(x1 .- x2)

                    if d < 2 * h
                        rho += cubic_spline_kernel(d, h)
                    end
                end
            end
        end

        densities[tid] = rho
    end

    makeindex!(; xss, props, cibuf, nvalid, grid, cell, backend, topo)
    kernel!(backend)(; ndrange=nvalid)

    @view densities[1:nvalid]
end

@inline function cubic_spline_kernel(r::Real, h::Real)
    q = r / h
    factor = 1 / (π * h^3)
    if q < 1
        factor * (1 - 1.5 * q^2 + 0.75 * q^3)
    elseif q < 2
        factor * 0.25 * (2 - q)^3
    else
        0.0
    end
end
