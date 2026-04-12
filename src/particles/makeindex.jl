xss_makeindex(::LCS.RKStage1, props::Properties) = props.xss
xss_makeindex(::LCS.RKStage2, props::Properties) = props.xss2

"""
    makeindex!(buf, state, config, topo, stage, iprofile)

Build cell-based spatial index for particles.
"""
function makeindex!(
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    stage::LCS.RKStage,
    iprofile::Integer,
)
    makeindex!(;
        xss=xss_makeindex(stage, buf.particles[iprofile].props),
        buf.particles[iprofile].props,
        cibuf=buf.particles[iprofile].ci,
        state.particles[iprofile].nvalid,
        config.grid,
        config.particles[iprofile].cell,
        config.backend,
        topo,
    )
end

function makeindex!(;
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
    (; perm, hashes, hasher, starts, stops, copy) = cibuf

    x0 = LCS.origins_l(topo)
    dxc = LCS.spacings(cell, grid)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x = xs[i], ys[i], zs[i]
            ci = Utils.unsafe_floor.(Int32, (x .- x0) ./ dxc) .+ 1
            hashes[i] = Utils.encode(hasher, ci)
        end
    end

    Utils.sortperm!(perm, hashes, nvalid; temp=reinterpret(eltype(perm), copy))
    Utils.sortbyperm!(hashes, perm, copy, nvalid)
    Utils.sortbyperm!(props, perm, copy, nvalid)

    h = Utils.eachhash(hasher)
    Parallel.foraxes(backend, (eachindex(h),)) do i
        @inbounds begin
            x = h[i]
            starts[x] = Utils.searchsortedfirst(hashes, x, nvalid)
            stops[x] = Utils.searchsortedlast(hashes, x, nvalid)
        end
    end

    LCS.EMPTY_LOG
end
