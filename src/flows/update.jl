Us_update(::LCS.RKStage1, fields::Fields) = fields.Us2
Us0_update(::LCS.RKStage1, fields::Fields) = fields.Us

Us_update(::LCS.RKStage2, fields::Fields) = fields.Us
Us0_update(::LCS.RKStage2, fields::Fields) = fields.Us

"""
    update!(buf, state, config, topo, stage)

Update velocity fields using explicit Runge-Kutta time integration.
"""
function update!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage)
    update!(;
        Us=Us_update(stage, buf.flow.fields),
        Us0=Us0_update(stage, buf.flow.fields),
        buf.flow.fields,
        dt=LCS.dtrk(stage, state),
        config.grid,
        config.backend,
        topo,
    )
end
function update!(;
    Us::Tuple3{Field},
    Us0::Tuple3{Field},
    fields::Fields,
    dt::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    U0, V0, W0 = Us0
    dUdt, dVdt, dWdt = fields.dUdts

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))
    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            du = $dUdt[i, j, k] * dt
            dv = $dVdt[i, j, k] * dt
            dw = $dWdt[i, j, k] * dt

            u = U0[i, j, k]
            v = V0[i, j, k]
            w = W0[i, j, k]

            U[i, j, k] = u + du
            V[i, j, k] = v + dv
            W[i, j, k] = w + dw
        end
    end

    LCS.EMPTY_LOG
end
