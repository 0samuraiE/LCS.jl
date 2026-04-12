function kwargs_init(::RandomFlow, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    (; buf.flow.fields.Us, buf.flow.fields, config.grid, topo)
end

function init!(init::RandomFlow; Us::Tuple3{Field}, fields::Fields, grid::LCS.Grid, topo::Topologies.Topology)
    U, V, W = Us
    P = fields.P

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))

    @offsetviews grid begin
        U_domain = U[domain...]
        V_domain = V[domain...]
        W_domain = W[domain...]
        P_domain = P[domain...]
    end

    s = init.seed + Topologies.linear_rank(topo)
    rng = Random.seed!(s)
    rand!(rng, U_domain)
    rand!(rng, V_domain)
    rand!(rng, W_domain)

    @. U_domain = (U_domain - 0.5) * 0.01
    @. V_domain = (V_domain - 0.5) * 0.01
    @. W_domain = (W_domain - 0.5) * 0.01
    @. P_domain = 0.0

    LCS.EMPTY_LOG
end
