@inline function hsmac_cellwise!(
    i::Integer,
    j::Integer,
    k::Integer,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Div2::Field,
    dti::Real,
    dxi::Real,
    dyi::Real,
    dzi::Real,
    dl::Real,
    omega::Real,
    grid::LCS.Grid,
)
    @inbounds @offsetviews grid begin
        dudx = dxi * LCS.diff1_o2_at_half(U[(i - 1):i, j, k])
        dvdy = dyi * LCS.diff1_o2_at_half(V[i, (j - 1):j, k])
        dwdz = dzi * LCS.diff1_o2_at_half(W[i, j, (k - 1):k])

        div = dudx + dvdy + dwdz
        dp = -omega * div * dl

        Div2[i, j, k] = div^2
        P[i, j, k] += dp * dti
        U[i - 1, j, k] -= dp * dxi
        U[i, j, k] += dp * dxi
        V[i, j - 1, k] -= dp * dyi
        V[i, j, k] += dp * dyi
        W[i, j, k - 1] -= dp * dzi
        W[i, j, k] += dp * dzi
    end
end

function kwargs_couple(
    ::HSMAC, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage
)
    (;
        Us=Us_couple(stage, buf.flow.fields),
        buf.flow.fields,
        buf.flow.scratches,
        halobuf=buf.halo,
        dt=LCS.dtrk(stage, state),
        config.grid,
        config.backend,
        topo,
    )
end

function couple!(
    couple::HSMAC;
    Us::Tuple3{Field},
    fields::Fields,
    scratches::Scratches,
    halobuf::Topologies.HaloBuffer,
    dt::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    P = fields.P
    Div2 = scratches[1]

    dxmin = minimum(LCS.spacings(grid))
    dti = 1 / dt

    if couple.time_blocking
        _couple_time_blocking!(couple, U, V, W, P, Div2, halobuf, dti, dxmin, grid, backend, topo)
    else
        _couple_naive!(couple, U, V, W, P, Div2, halobuf, dti, dxmin, grid, backend, topo)
    end
end

function _couple_naive!(
    couple::HSMAC,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Div2::Field,
    halobuf::Topologies.HaloBuffer,
    dti::Real,
    dxmin::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)
    dl = 1 / (2 * (dxi^2 + dyi^2 + dzi^2))
    omega = couple.omega

    itr, div = 0, 0.0
    while itr < couple.itrmax
        domain = map(n -> 1:(n + 1), LCS.dims_l(grid, topo))
        if couple.overlap
            core = map(n -> 2:n, LCS.dims_l(grid, topo))

            Topologies.synchalo!(Topologies.FullHalo(1), (U, V, W), halobuf, grid, backend, topo) do
                Parallel.foraxes(backend, core; couple.coloring) do i, j, k
                    hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                end
            end

            for face in Topologies.faces(domain, core)
                Parallel.foraxes(backend, face; couple.coloring) do i, j, k
                    hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                end
            end
        else
            Topologies.synchalo!(Topologies.FullHalo(1), (U, V, W, P), halobuf, grid, backend, topo)

            Parallel.foraxes(backend, domain; couple.coloring) do i, j, k
                hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
            end
        end

        div = @offsetviews grid sqrt(Topologies.allmean(Div2[domain...], topo))

        if div * dxmin < couple.epsp && couple.itrmin <= itr
            break
        end
        itr += 1
    end

    (; itr, div)
end

function _couple_time_blocking!(
    couple::HSMAC,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Div2::Field,
    halobuf::Topologies.HaloBuffer,
    dti::Real,
    dxmin::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)
    dl = 1 / (2 * (dxi^2 + dyi^2 + dzi^2))
    omega = couple.omega

    itr, div = 0, 0.0
    while itr < couple.itrmax
        domain1 = map(n -> -1:(n + 3), LCS.dims_l(grid, topo))
        domain2 = map(n -> 0:(n + 2), LCS.dims_l(grid, topo))
        domain3 = map(n -> 1:(n + 1), LCS.dims_l(grid, topo))

        if couple.overlap
            core = map(n -> 2:n, LCS.dims_l(grid, topo))

            for domain in (domain1, domain2, domain3)
                if domain === domain1
                    Topologies.synchalo!(Topologies.FullHalo(grid), (U, V, W), halobuf, grid, backend, topo) do
                        Parallel.foraxes(backend, core; couple.coloring) do i, j, k
                            hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                        end
                    end

                    for face in Topologies.faces(domain1, core)
                        Parallel.foraxes(backend, face; couple.coloring) do i, j, k
                            hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                        end
                    end
                else
                    Parallel.foraxes(backend, domain; couple.coloring) do i, j, k
                        hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                    end
                end

                div = @offsetviews grid sqrt(Topologies.allmean(Div2[domain3...], topo))
                if div * dxmin < couple.epsp && couple.itrmin <= itr
                    return (; itr, div)
                end
                itr += 1
            end
        else
            Topologies.synchalo!(Topologies.FullHalo(grid), (U, V, W), halobuf, grid, backend, topo)

            for domain in (domain1, domain2, domain3)
                Parallel.foraxes(backend, domain; couple.coloring) do i, j, k
                    hsmac_cellwise!(i, j, k, U, V, W, P, Div2, dti, dxi, dyi, dzi, dl, omega, grid)
                end

                div = @offsetviews grid sqrt(Topologies.allmean(Div2[domain3...], topo))
                if div * dxmin < couple.epsp && couple.itrmin <= itr
                    return (; itr, div)
                end
                itr += 1
            end
        end
    end

    (; itr, div)
end
