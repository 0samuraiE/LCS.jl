Us_rhs(::LCS.RKStage1, fields::Fields) = fields.Us
Us_rhs(::LCS.RKStage2, fields::Fields) = fields.Us2

"""
    rhs!(buf, state, config, topo, stage)

Compute right-hand side of incompressible Navier-Stokes equations.
"""
function rhs!(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage)
    rhs!(;
        Us=Us_rhs(stage, buf.flow.fields),
        buf.flow.fields,
        buf.flow.scratches,
        config.grid,
        config.flow.params,
        config.backend,
        topo,
    )
end
function rhs!(;
    Us::Tuple3{Field},
    fields::Fields,
    scratches::Scratches,
    grid::LCS.Grid,
    params::FlowParams,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    P = fields.P
    dUdt, dVdt, dWdt = fields.dUdts
    Jx1, Jy1, Jz1, Jx3, Jy3, Jz3 = scratches

    domain1 = map(n -> -1:(n + 1), LCS.dims_l(grid, topo))
    domain2 = map(n -> 1:n, LCS.dims_l(grid, topo))

    rhs_x!(dUdt, U, V, W, P, Jx1, Jy1, Jz1, Jx3, Jy3, Jz3, grid, params, backend, domain1, domain2)
    rhs_y!(dVdt, U, V, W, P, Jx1, Jy1, Jz1, Jx3, Jy3, Jz3, grid, params, backend, domain1, domain2)
    rhs_z!(dWdt, U, V, W, P, Jx1, Jy1, Jz1, Jx3, Jy3, Jz3, grid, params, backend, domain1, domain2)

    LCS.EMPTY_LOG
end

function rhs_x!(
    dUdt::Field,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Jx1::Field,
    Jy1::Field,
    Jz1::Field,
    Jx3::Field,
    Jy3::Field,
    Jz3::Field,
    grid::LCS.Grid,
    params::FlowParams,
    backend::KA.Backend,
    domain1::Tuple3{AbstractUnitRange},
    domain2::Tuple3{AbstractUnitRange},
)
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)
    Rei = 1 / Flows.Re(params)

    Parallel.foraxes(backend, domain1) do i, j, k
        @inbounds @offsetviews grid begin
            uc = LCS.interp_o4_at_half(U[(i - 1):(i + 2), j, k])
            vc = LCS.interp_o4_at_half(V[(i - 1):(i + 2), j, k])
            wc = LCS.interp_o4_at_half(W[(i - 1):(i + 2), j, k])

            ux1 = LCS.interp_o2_at_half(U[i:(i + 1), j, k])
            uy1 = LCS.interp_o2_at_half(U[i, j:(j + 1), k])
            uz1 = LCS.interp_o2_at_half(U[i, j, k:(k + 1)])

            ux3 = LCS.interp_o2_at_half(U[(i - 1):3:(i + 2), j, k])
            uy3 = LCS.interp_o2_at_half(U[i, (j - 1):3:(j + 2), k])
            uz3 = LCS.interp_o2_at_half(U[i, j, (k - 1):3:(k + 2)])

            if CartesianIndex(i, j, k) in CartesianIndices(domain2)
                diff =
                    (
                        dxi^2 * LCS.diff2_o4_at_node(U[(i - 2):(i + 2), j, k]) +
                        dyi^2 * LCS.diff2_o4_at_node(U[i, (j - 2):(j + 2), k]) +
                        dzi^2 * LCS.diff2_o4_at_node(U[i, j, (k - 2):(k + 2)])
                    ) * Rei

                pres = dxi * LCS.diff1_o2_at_half(P[i:(i + 1), j, k])

                $dUdt[i, j, k] = diff - pres
            end

            Jx1[i, j, k] = uc * ux1
            Jy1[i, j, k] = vc * uy1
            Jz1[i, j, k] = wc * uz1

            Jx3[i, j, k] = uc * ux3
            Jy3[i, j, k] = vc * uy3
            Jz3[i, j, k] = wc * uz3
        end
    end

    Parallel.foraxes(backend, domain2) do i, j, k
        @inbounds @offsetviews grid begin
            advect =
                dxi * LCS.sum_flux_o4(Jx1[(i - 1):i, j, k], Jx3[(i - 2):3:(i + 1), j, k]) +
                dyi * LCS.sum_flux_o4(Jy1[i, (j - 1):j, k], Jy3[i, (j - 2):3:(j + 1), k]) +
                dzi * LCS.sum_flux_o4(Jz1[i, j, (k - 1):k], Jz3[i, j, (k - 2):3:(k + 1)])

            $dUdt[i, j, k] -= advect
        end
    end
end

function rhs_y!(
    dVdt::Field,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Jx1::Field,
    Jy1::Field,
    Jz1::Field,
    Jx3::Field,
    Jy3::Field,
    Jz3::Field,
    grid::LCS.Grid,
    params::FlowParams,
    backend::KA.Backend,
    domain1::Tuple3{AbstractUnitRange},
    domain2::Tuple3{AbstractUnitRange},
)
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)
    Rei = 1 / Flows.Re(params)

    Parallel.foraxes(backend, domain1) do i, j, k
        @inbounds @offsetviews grid begin
            uc = LCS.interp_o4_at_half(U[i, (j - 1):(j + 2), k])
            vc = LCS.interp_o4_at_half(V[i, (j - 1):(j + 2), k])
            wc = LCS.interp_o4_at_half(W[i, (j - 1):(j + 2), k])

            vx1 = LCS.interp_o2_at_half(V[i:(i + 1), j, k])
            vy1 = LCS.interp_o2_at_half(V[i, j:(j + 1), k])
            vz1 = LCS.interp_o2_at_half(V[i, j, k:(k + 1)])

            vx3 = LCS.interp_o2_at_half(V[(i - 1):3:(i + 2), j, k])
            vy3 = LCS.interp_o2_at_half(V[i, (j - 1):3:(j + 2), k])
            vz3 = LCS.interp_o2_at_half(V[i, j, (k - 1):3:(k + 2)])

            if CartesianIndex(i, j, k) in CartesianIndices(domain2)
                diff =
                    (
                        dxi^2 * LCS.diff2_o4_at_node(V[(i - 2):(i + 2), j, k]) +
                        dyi^2 * LCS.diff2_o4_at_node(V[i, (j - 2):(j + 2), k]) +
                        dzi^2 * LCS.diff2_o4_at_node(V[i, j, (k - 2):(k + 2)])
                    ) * Rei

                pres = dyi * LCS.diff1_o2_at_half(P[i, j:(j + 1), k])

                $dVdt[i, j, k] = diff - pres
            end

            Jx1[i, j, k] = uc * vx1
            Jy1[i, j, k] = vc * vy1
            Jz1[i, j, k] = wc * vz1

            Jx3[i, j, k] = uc * vx3
            Jy3[i, j, k] = vc * vy3
            Jz3[i, j, k] = wc * vz3
        end
    end

    Parallel.foraxes(backend, domain2) do i, j, k
        @inbounds @offsetviews grid begin
            advect =
                dxi * LCS.sum_flux_o4(Jx1[(i - 1):i, j, k], Jx3[(i - 2):3:(i + 1), j, k]) +
                dyi * LCS.sum_flux_o4(Jy1[i, (j - 1):j, k], Jy3[i, (j - 2):3:(j + 1), k]) +
                dzi * LCS.sum_flux_o4(Jz1[i, j, (k - 1):k], Jz3[i, j, (k - 2):3:(k + 1)])

            $dVdt[i, j, k] -= advect
        end
    end
end

function rhs_z!(
    dWdt::Field,
    U::Field,
    V::Field,
    W::Field,
    P::Field,
    Jx1::Field,
    Jy1::Field,
    Jz1::Field,
    Jx3::Field,
    Jy3::Field,
    Jz3::Field,
    grid::LCS.Grid,
    params::FlowParams,
    backend::KA.Backend,
    domain1::Tuple3{AbstractUnitRange},
    domain2::Tuple3{AbstractUnitRange},
)
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)
    Rei = 1 / Flows.Re(params)

    Parallel.foraxes(backend, domain1) do i, j, k
        @inbounds @offsetviews grid begin
            uc = LCS.interp_o4_at_half(U[i, j, (k - 1):(k + 2)])
            vc = LCS.interp_o4_at_half(V[i, j, (k - 1):(k + 2)])
            wc = LCS.interp_o4_at_half(W[i, j, (k - 1):(k + 2)])

            wx1 = LCS.interp_o2_at_half(W[i:(i + 1), j, k])
            wy1 = LCS.interp_o2_at_half(W[i, j:(j + 1), k])
            wz1 = LCS.interp_o2_at_half(W[i, j, k:(k + 1)])

            wx3 = LCS.interp_o2_at_half(W[(i - 1):3:(i + 2), j, k])
            wy3 = LCS.interp_o2_at_half(W[i, (j - 1):3:(j + 2), k])
            wz3 = LCS.interp_o2_at_half(W[i, j, (k - 1):3:(k + 2)])

            if CartesianIndex(i, j, k) in CartesianIndices(domain2)
                diff =
                    (
                        dxi^2 * LCS.diff2_o4_at_node(W[(i - 2):(i + 2), j, k]) +
                        dyi^2 * LCS.diff2_o4_at_node(W[i, (j - 2):(j + 2), k]) +
                        dzi^2 * LCS.diff2_o4_at_node(W[i, j, (k - 2):(k + 2)])
                    ) * Rei

                pres = dzi * LCS.diff1_o2_at_half(P[i, j, k:(k + 1)])

                $dWdt[i, j, k] = diff - pres
            end

            Jx1[i, j, k] = uc * wx1
            Jy1[i, j, k] = vc * wy1
            Jz1[i, j, k] = wc * wz1

            Jx3[i, j, k] = uc * wx3
            Jy3[i, j, k] = vc * wy3
            Jz3[i, j, k] = wc * wz3
        end
    end

    Parallel.foraxes(backend, domain2) do i, j, k
        @inbounds @offsetviews grid begin
            advect =
                dxi * LCS.sum_flux_o4(Jx1[(i - 1):i, j, k], Jx3[(i - 2):3:(i + 1), j, k]) +
                dyi * LCS.sum_flux_o4(Jy1[i, (j - 1):j, k], Jy3[i, (j - 2):3:(j + 1), k]) +
                dzi * LCS.sum_flux_o4(Jz1[i, j, (k - 1):k], Jz3[i, j, (k - 2):3:(k + 1)])

            $dWdt[i, j, k] -= advect
        end
    end
end
