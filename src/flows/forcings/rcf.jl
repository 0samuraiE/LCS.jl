function forcing!(
    forcing::EnergyPreserveRCF,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    stage::LCS.RKStage2,
)
    forcing!(
        forcing;
        buf.flow.fields.Us,
        Power=buf.flow.scratches[1],
        forcingbuf=buf.flow.forcing,
        dt=LCS.dtrk(stage, state),
        config.grid,
        config.backend,
        topo,
    )
end

function forcing!(
    forcing::EnergyPreserveRCF;
    Us::Tuple3{Field},
    Power::Field,
    forcingbuf::RCFBuffer,
    dt::Real,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us

    U_hat, V_hat, W_hat = forcingbuf.Us_hat
    U_coarse_g, V_coarse_g, W_coarse_g = forcingbuf.Us_coarse_g

    dti = 1 / dt

    phi = Flows.phase_shift(forcing, grid)
    e10 = forcing.energy
    kmin = forcing.kmin
    kmax = forcing.kmax

    boxmean!(Topologies.processing(topo); Us, forcingbuf, grid, forcing, backend, topo)

    Utils.rfft!(U_hat, forcingbuf.plan, U_coarse_g)
    Utils.rfft!(V_hat, forcingbuf.plan, V_coarse_g)
    Utils.rfft!(W_hat, forcingbuf.plan, W_coarse_g)

    for k0 in kmin:kmax
        en = FT.enespe(U_hat, V_hat, W_hat, k0)
        en0 = enespe0(k0, e10)
        r = sqrt(en0 / max(en, 1E-8)) - 1.0
        for kz in (-kmax):kmax, ky in (-kmax):kmax, kx in 0:kmax
            r2 = kx^2 + ky^2 + kz^2

            if (k0 - 0.5)^2 <= r2 < (k0 + 0.5)^2
                i = FT.fftindex.((kx, ky, kz), size(U_hat))
                U_hat[i...] *= r
                V_hat[i...] *= r
                W_hat[i...] *= r
            end
        end
    end

    U_hat_cropped = FT.crop(U_hat, kmax; forcingbuf.temp)
    V_hat_cropped = FT.crop(V_hat, kmax; forcingbuf.temp)
    W_hat_cropped = FT.crop(W_hat, kmax; forcingbuf.temp)

    vkmin, vkmax = Val(kmin), Val(kmax)
    nxi_g, nyi_g, nzi_g = 1 ./ LCS.dims_g(grid)
    dims_l = LCS.dims_l(grid, topo)
    domain = map(n -> 1:n, dims_l)

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsets grid begin
            i_g, j_g, k_g = LCS.index_g((i, j, k), dims_l, topo)

            exp1k = cispi(2 * (k_g - phi) * nzi_g)
            exp1j = cispi(2 * (j_g - phi) * nyi_g)
            exp1i = cispi(2 * (i_g - phi) * nxi_g)

            usrc = FT.sirdft(U_hat_cropped, exp1i, exp1j, exp1k, Parallel.unwrap(vkmin), Parallel.unwrap(vkmax))
            vsrc = FT.sirdft(V_hat_cropped, exp1i, exp1j, exp1k, Parallel.unwrap(vkmin), Parallel.unwrap(vkmax))
            wsrc = FT.sirdft(W_hat_cropped, exp1i, exp1j, exp1k, Parallel.unwrap(vkmin), Parallel.unwrap(vkmax))

            local u = U[i, j, k]
            local v = V[i, j, k]
            local w = W[i, j, k]
            local p = (u * usrc + v * vsrc + w * wsrc) * dti

            U[i, j, k] += usrc
            V[i, j, k] += vsrc
            W[i, j, k] += wsrc
            Power[i, j, k] = p
        end
    end

    p = @offsetviews grid Topologies.allmean(Power[domain...], topo)
    (; p, e10)
end

function boxmean!(
    ::Topologies.SingleProcessing;
    Us::Tuple3{Field},
    forcingbuf::RCFBuffer,
    grid::LCS.Grid,
    forcing::RCForcing,
    kwargs...,
)
    (; Us_coarse_l, Us_coarse_g) = forcingbuf

    halo_size = grid.halo_size
    filter_dims = Flows.filter_dims(forcing, grid)

    U, V, W = Us
    U_coarse_g, V_coarse_g, W_coarse_g = Us_coarse_g
    U_coarse_l = view(Us_coarse_l, :, :, :, 1)
    V_coarse_l = view(Us_coarse_l, :, :, :, 2)
    W_coarse_l = view(Us_coarse_l, :, :, :, 3)

    boxmean_l!(U_coarse_l, U, halo_size, filter_dims)
    boxmean_l!(V_coarse_l, V, halo_size, filter_dims)
    boxmean_l!(W_coarse_l, W, halo_size, filter_dims)

    copyto!(U_coarse_g, U_coarse_l)
    copyto!(V_coarse_g, V_coarse_l)
    copyto!(W_coarse_g, W_coarse_l)
end

function boxmean!(
    ::Topologies.MultiProcessing;
    Us::Tuple3{Field},
    forcingbuf::RCFBuffer,
    grid::LCS.Grid,
    forcing::RCForcing,
    backend::KA.Backend,
    topo::Topologies.Topology,
    kwargs...,
)
    (; Us_coarse_l, Us_coarse_g, sendbuf, recvbuf) = forcingbuf
    !isnothing(sendbuf) && !isnothing(recvbuf) ||
        throw(ArgumentError("sendbuf and recvbuf must be provided for multiprocessing"))

    halo_size = grid.halo_size
    filter_dims = Flows.filter_dims(forcing, grid)
    coarse_dims_l = Flows.coarse_dims_l(forcing, topo)

    comm = Topologies.comm()
    cart_rank = Topologies.cart_rank(topo)
    proc_size = Topologies.proc_size(topo)

    U, V, W = Us
    U_coarse_g, V_coarse_g, W_coarse_g = Us_coarse_g
    U_coarse_l = view(Us_coarse_l, :, :, :, 1)
    V_coarse_l = view(Us_coarse_l, :, :, :, 2)
    W_coarse_l = view(Us_coarse_l, :, :, :, 3)

    boxmean_l!(U_coarse_l, U, halo_size, filter_dims)
    boxmean_l!(V_coarse_l, V, halo_size, filter_dims)
    boxmean_l!(W_coarse_l, W, halo_size, filter_dims)

    copyto!(sendbuf, Us_coarse_l)
    KA.synchronize(backend)
    MPI.Allgather!(sendbuf, recvbuf, comm)

    recvbuf = Utils.reshape(recvbuf, size(Us_coarse_l)..., proc_size)
    @views for cart_rank in Topologies.each_cart_rank(topo)
        indices = LCS.indices_l(coarse_dims_l, cart_rank)
        rank = Topologies.cart_to_linear_rank(topo, cart_rank)

        copyto!(U_coarse_g[indices...], recvbuf[:, :, :, 1, rank + 1])
        copyto!(V_coarse_g[indices...], recvbuf[:, :, :, 2, rank + 1])
        copyto!(W_coarse_g[indices...], recvbuf[:, :, :, 3, rank + 1])
    end
end

function boxmean_l!(U_coarse_l::Field, U::Field, halo_size::Integer, filter_dims::Integer3)
    U = Utils.reshape(U, size(U)..., 1, 1)
    U_coarse_l = Utils.reshape(U_coarse_l, size(U_coarse_l)..., 1, 1)

    k = stride = filter_dims
    padding = -halo_size
    meanpool!(U_coarse_l, U, PoolDims(size(U), k; stride, padding))
end

function enespe0(k0::Integer, e10::Real)
    c1 = FT.kcount(1)
    ck = FT.kcount(k0)
    e10 * k0^(-5 / 3) * (1.5^3 - 0.5^3) / ((k0 + 0.5)^3 - (k0 - 0.5)^3) * ck / c1
end
