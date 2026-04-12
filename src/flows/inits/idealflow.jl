function kwargs_init(::IdealFlow, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    (; buf.flow.fields.Us, buf.flow.fields, config.grid, config.backend, topo)
end

function init!(
    ::IdealFlow; Us::Tuple3{Field}, fields::Fields, grid::LCS.Grid, backend::KA.Backend, topo::Topologies.Topology
)
    U, V, W = Us
    P = fields.P

    dims_g = LCS.dims_g(grid)
    allequal(dims_g) || throw(ArgumentError("grid dimensions must be cubic, got $dims_g"))
    n = first(dims_g)

    dx, dy, dz = LCS.spacings(grid)
    dims_l = LCS.dims_l(grid, topo)
    domain = map(n -> 1:n, dims_l)

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            i_g, j_g, k_g = LCS.index_g((i, j, k), dims_l, topo)

            x_g = dx * (i_g - 0.5)
            y_g = dy * (j_g - 0.5)
            z_g = dz * (k_g - 0.5)

            al, bl, cl = 1.0, 1.0, 0.0
            Al, Bl, Cl = sqrt(1.5) * 2.0, -sqrt(1.5) * 2.0, 0.0

            as, bs, cs = 1.0 * (n / 2), 1.0 * (n / 2), 1.0 * (n / 2)
            As, Bs, Cs = sqrt(1.5e-2) * 2.0, 0.0, -sqrt(1.5e-2) * 2.0

            ul = Al * cos(al * (x_g + 0.5 * dx)) * sin(bl * y_g)
            vl = Bl * sin(al * x_g) * cos(bl * (y_g + 0.5 * dy))
            wl = Cl * sin(al * x_g) * sin(bl * y_g)

            us = As * cos(as * (x_g + 0.5 * dx)) * sin(bs * y_g) * sin(cs * z_g)
            vs = Bs * sin(as * x_g) * cos(bs * (y_g + 0.5 * dy)) * sin(cs * z_g)
            ws = Cs * sin(as * x_g) * sin(bs * y_g) * cos(cs * (z_g + 0.5 * dz))

            U[i, j, k] = ul + us
            V[i, j, k] = vl + vs
            W[i, j, k] = wl + ws
            P[i, j, k] = 0.0
        end
    end

    LCS.EMPTY_LOG
end
