function load_coarse_flow_fields!(
    file::String, buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology
)
    U, V, W = buf.flow.fields.Us
    P = buf.flow.fields.P
    A, B, C, D = buf.flow.scratches

    dims_l = LCS.dims_l(config.grid, topo)
    dims_l_coarse = div.(dims_l, 2)
    halo_size = config.grid.halo_size
    backend = config.backend

    subindices = map(n -> 1:(n + 2 * halo_size), dims_l_coarse)
    @views begin
        Uc = A[subindices...]
        Vc = B[subindices...]
        Wc = C[subindices...]
        Pc = D[subindices...]
    end

    _h5open(file, "r", topo) do h
        g = h[TAG_FLOW]
        indices = map(n -> 1:n, dims_l_coarse)

        @offsetviews config.grid begin
            @log backend "io/read/flow/coarse/U" read_field(g, TAG_FLOW_U, Uc[indices...], topo)
            @log backend "io/read/flow/coarse/V" read_field(g, TAG_FLOW_V, Vc[indices...], topo)
            @log backend "io/read/flow/coarse/W" read_field(g, TAG_FLOW_W, Wc[indices...], topo)
            @log backend "io/read/flow/coarse/P" read_field(g, TAG_FLOW_P, Pc[indices...], topo)
        end
    end

    KA.synchronize(backend)

    grid_coarse = LCS.Grid(div.(config.grid.dims, 2), halo_size)
    Topologies.synchalo!(Topologies.FullHalo(1), (Uc, Vc, Wc, Pc), buf.halo, grid_coarse, backend, topo)

    trilinear_upsample2x!(U, Uc, dims_l_coarse, config.grid, backend)
    trilinear_upsample2x!(V, Vc, dims_l_coarse, config.grid, backend)
    trilinear_upsample2x!(W, Wc, dims_l_coarse, config.grid, backend)
    trilinear_upsample2x!(P, Pc, dims_l_coarse, config.grid, backend)
end

function trilinear_upsample2x!(
    A_fine::AbstractArray, A_coarse::AbstractArray, dims_l_coarse::Integer3, grid::LCS.Grid, backend::KA.Backend
)
    domain = map(n -> 1:n, dims_l_coarse)
    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            ii, jj, kk = 2i - 1, 2j - 1, 2k - 1
            A_fine[ii, jj, kk] = A_coarse[i, j, k]

            ii, jj, kk = 2i, 2j - 1, 2k - 1
            A_fine[ii, jj, kk] = (1 / 2) * (A_coarse[i, j, k] + A_coarse[i + 1, j, k])

            ii, jj, kk = 2i - 1, 2j, 2k - 1
            A_fine[ii, jj, kk] = (1 / 2) * (A_coarse[i, j, k] + A_coarse[i, j + 1, k])

            ii, jj, kk = 2i - 1, 2j - 1, 2k
            A_fine[ii, jj, kk] = (1 / 2) * (A_coarse[i, j, k] + A_coarse[i, j, k + 1])

            ii, jj, kk = 2i, 2j, 2k - 1
            A_fine[ii, jj, kk] =
                (1 / 4) *
                (A_coarse[i, j, k] + A_coarse[i + 1, j, k] + A_coarse[i, j + 1, k] + A_coarse[i + 1, j + 1, k])

            ii, jj, kk = 2i, 2j - 1, 2k
            A_fine[ii, jj, kk] =
                (1 / 4) *
                (A_coarse[i, j, k] + A_coarse[i + 1, j, k] + A_coarse[i, j, k + 1] + A_coarse[i + 1, j, k + 1])

            ii, jj, kk = 2i - 1, 2j, 2k
            A_fine[ii, jj, kk] =
                (1 / 4) *
                (A_coarse[i, j, k] + A_coarse[i, j + 1, k] + A_coarse[i, j, k + 1] + A_coarse[i, j + 1, k + 1])

            ii, jj, kk = 2i, 2j, 2k
            A_fine[ii, jj, kk] =
                (1 / 8) * (
                    A_coarse[i, j, k] +
                    A_coarse[i + 1, j, k] +
                    A_coarse[i, j + 1, k] +
                    A_coarse[i, j, k + 1] +
                    A_coarse[i + 1, j + 1, k] +
                    A_coarse[i + 1, j, k + 1] +
                    A_coarse[i, j + 1, k + 1] +
                    A_coarse[i + 1, j + 1, k + 1]
                )
        end
    end
end
