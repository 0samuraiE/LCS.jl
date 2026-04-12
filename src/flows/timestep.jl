"""
    timestep(buf, state, config, topo)

Compute adaptive time step size for stable flow simulation.
"""
function timestep(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    timestep(;
        buf.flow.fields.Us,
        buf.flow.scratches,
        halobuf=buf.halo,
        config.grid,
        ts=config.timestep,
        config.flow.params,
        config.backend,
        topo,
    )
end
function timestep(;
    Us::Tuple3{Field},
    scratches::Scratches,
    halobuf::Topologies.HaloBuffer,
    grid::LCS.Grid,
    ts::LCS.TimeStep,
    params::FlowParams,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    U_interp, V_interp, W_interp = scratches

    Topologies.synchalo!(Topologies.FaceHalo(grid), Us, halobuf, grid, backend, topo)

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            uc = LCS.interp_o4_at_half(U[(i - 2):(i + 1), j, k])
            vc = LCS.interp_o4_at_half(V[i, (j - 2):(j + 1), k])
            wc = LCS.interp_o4_at_half(W[i, j, (k - 2):(k + 1)])

            U_interp[i, j, k] = uc
            V_interp[i, j, k] = vc
            W_interp[i, j, k] = wc
        end
    end

    @inbounds @offsetviews grid begin
        umax = maximum(abs, U_interp[domain...]) + 1e-10
        vmax = maximum(abs, V_interp[domain...]) + 1e-10
        wmax = maximum(abs, W_interp[domain...]) + 1e-10
    end

    dx, dy, dz = LCS.spacings(grid)
    cfl = ts.cfl
    Re = Flows.Re(params)
    cands = (#
        ts.dtmax,
        cfl * dx / umax,
        cfl * dy / vmax,
        cfl * dz / wmax,
        0.5 * (dx / 2)^2 / (1 / Re),
    )

    Topologies.allminimum(cands, topo)
end

"""
    timestep!(buf, state, config, topo)

Update time step in state based on stability constraints.
"""
function timestep!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    state.dt = timestep(buf, state, config, topo)
end
