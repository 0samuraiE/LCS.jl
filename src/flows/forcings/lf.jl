"""
    forcing!(forcing, buf, state, config, topo, stage)

Apply ConstantPowerLF forcing: scales velocity fluctuations to inject constant power.
"""
function forcing!(
    forcing::ConstantPowerLF,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    stage::LCS.RKStage2,
)
    forcing!(
        forcing;
        buf.flow.fields.Us,
        U2=buf.flow.scratches[1],
        halobuf=buf.halo,
        state.t,
        dt=LCS.dtrk(stage, state),
        config.grid,
        config.backend,
        topo,
    )
end

function forcing!(
    forcing::ConstantPowerLF;
    Us::Tuple3{Field},
    U2::Field,
    halobuf::Topologies.HaloBuffer,
    t::Real,
    dt::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    p0 = forcing.power

    Topologies.synchalo!(Topologies.FullHalo(grid), Us, halobuf, grid, backend, topo)

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            local u = LCS.interp_o2_at_half(U[(i - 1):i, j, k])
            local v = LCS.interp_o2_at_half(V[i, (j - 1):j, k])
            local w = LCS.interp_o2_at_half(W[i, j, (k - 1):k])

            local u2 = u^2 + v^2 + w^2
            U2[i, j, k] = u2
        end
    end

    @offsetviews grid begin
        u2 = Topologies.allmean(U2[domain...], topo)
        u  = Topologies.allmean(U[domain...], topo)
        v  = Topologies.allmean(V[domain...], topo)
        w  = Topologies.allmean(W[domain...], topo)
    end

    e = 0.5 * (u2 - u^2 - v^2 - w^2)
    A = p0 / (2 * e)

    Power = U2
    p = _linear_forcing!(Us, Power, A, u, v, w, dt, grid, backend, topo)

    (; p)
end

function _linear_forcing!(
    Us::Tuple3{Field},
    Power::Field,
    A::Real,
    u_mean::Real,
    v_mean::Real,
    w_mean::Real,
    dt::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))
    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            usrc = A * (U[i, j, k] - u_mean)
            vsrc = A * (V[i, j, k] - v_mean)
            wsrc = A * (W[i, j, k] - w_mean)

            u = U[i, j, k]
            v = V[i, j, k]
            w = W[i, j, k]
            p = u * usrc + v * vsrc + w * wsrc

            U[i, j, k] = u + usrc * dt
            V[i, j, k] = v + vsrc * dt
            W[i, j, k] = w + wsrc * dt
            Power[i, j, k] = p
        end
    end

    @offsetviews grid Topologies.allmean(Power[domain...], topo)
end
