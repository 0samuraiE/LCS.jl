struct ExchangeSummary
    mean :: Int
    min  :: Int
    max  :: Int
end
@composite ExchangeSummary

function PPrint.pprint(io::IO, s::ExchangeSummary)
    PPrint.pprint(io, s.mean)
    print(io, "  ")
    PPrint.pprint(io, s.max)
    print(io, "  ")
    PPrint.pprint(io, s.min)
end

function ExchangeSummary(x::Integer, topo::Topologies.Topology)
    ExchangeSummary(
        ceil(Int, Topologies.allmean(x, topo)), Topologies.allminimum(x, topo), Topologies.allmaximum(x, topo)
    )
end

xss_boundary(::LCS.RKStage1, props::Properties) = props.xss
xss_boundary(::LCS.RKStage2, props::Properties) = props.xss2

xss2_boundary(::LCS.RKStage1, props::Properties) = props.xss2
xss2_boundary(::LCS.RKStage2, props::Properties) = props.xss

"""
    boundary!(buf, state, config, topo, stage, iprofile)

Apply boundary conditions and exchange particles between MPI processes.
"""
function boundary!(
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    stage::LCS.RKStage,
    iprofile::Integer,
)
    (; nvalid, nsend, nrecv) = _boundary!(;
        xss=xss_boundary(stage, buf.particles[iprofile].props),
        xss2=xss2_boundary(stage, buf.particles[iprofile].props),
        buf.particles[iprofile].props,
        commbuf=buf.particles[iprofile].comm,
        state.particles[iprofile].nvalid,
        config.particles[iprofile].population,
        config.backend,
        topo,
    )

    state.particles[iprofile].nvalid = nvalid

    (;#
        nvalid=ExchangeSummary(nvalid, topo),
        nrecvs=ExchangeSummary(nrecv, topo),
        nsends=ExchangeSummary(nsend, topo),
    )
end

function _boundary!(;
    xss::Tuple3{Property{<:Real}},
    xss2::Tuple3{Property{<:Real}},
    props::Properties,
    commbuf::Union{CommBuffer,Nothing},
    nvalid::Integer,
    population::Population,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    nsend, nrecv = 0, 0
    if Topologies.is_multi_processing(topo)
        ntotal = Particles.ntotal_l(population, topo)
        (; nvalid, nsend, nrecv) = exchange_props!(; xss, props, commbuf, ntotal, nvalid, backend, topo)
    end
    apply_periodic_boundary!(; xss, xss2, nvalid, backend)

    (; nvalid, nsend, nrecv)
end

function apply_periodic_boundary!(;
    xss::Tuple3{Property{<:Real}}, xss2::Tuple3{Property{<:Real}}, nvalid::Integer, backend::KA.Backend
)
    xs, ys, zs = xss
    xs2, ys2, zs2 = xss2

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x, y, z = xs[i], ys[i], zs[i]
            x2, y2, z2 = xs2[i], ys2[i], zs2[i]

            if x < LCS.DOMAIN_ORIGIN
                xs[i] = x + LCS.DOMAIN_LENGTH
                xs2[i] = x2 + LCS.DOMAIN_LENGTH
            end
            if x >= LCS.DOMAIN_ORIGIN + LCS.DOMAIN_LENGTH
                xs[i] = x - LCS.DOMAIN_LENGTH
                xs2[i] = x2 - LCS.DOMAIN_LENGTH
            end

            if y < LCS.DOMAIN_ORIGIN
                ys[i] = y + LCS.DOMAIN_LENGTH
                ys2[i] = y2 + LCS.DOMAIN_LENGTH
            end
            if y >= LCS.DOMAIN_ORIGIN + LCS.DOMAIN_LENGTH
                ys[i] = y - LCS.DOMAIN_LENGTH
                ys2[i] = y2 - LCS.DOMAIN_LENGTH
            end

            if z < LCS.DOMAIN_ORIGIN
                zs[i] = z + LCS.DOMAIN_LENGTH
                zs2[i] = z2 + LCS.DOMAIN_LENGTH
            end
            if z >= LCS.DOMAIN_ORIGIN + LCS.DOMAIN_LENGTH
                zs[i] = z - LCS.DOMAIN_LENGTH
                zs2[i] = z2 - LCS.DOMAIN_LENGTH
            end
        end
    end
end

function exchange_props!(;
    xss::Tuple3{Property{<:Real}},
    props::Properties,
    commbuf::CommBuffer,
    ntotal::Integer,
    nvalid::Integer,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    origins_l = LCS.origins_l(topo)
    ends_l = LCS.ends_l(topo)

    rank_self = Topologies.linear_rank(topo)
    neighbors = topo.neighbor_linear_ranks
    rank_rights = neighbors[1], neighbors[3], neighbors[5]
    rank_lefts = neighbors[2], neighbors[4], neighbors[6]

    nsend_total, nrecv_total = 0, 0
    for (axis, (xs, x_left, x_right, rank_left, rank_right)) in
        enumerate(zip(xss, origins_l, ends_l, rank_lefts, rank_rights))
        if rank_self == rank_left && rank_self == rank_right
            continue
        end

        (; nvalid, nsend, nrecv) = exchange_props_along_xaxis!(;
            xs, props, commbuf, backend, axis, rank_left, rank_right, nvalid, ntotal, x_left, x_right
        )
        nsend_total += nsend
        nrecv_total += nrecv
    end

    (; nvalid, nsend=nsend_total, nrecv=nrecv_total)
end

function exchange_props_along_xaxis!(;
    xs::Property{<:Real},
    props::Properties,
    commbuf::CommBuffer,
    backend::KA.Backend,
    axis::Integer,
    rank_left::Integer,
    rank_right::Integer,
    nvalid::Integer,
    ntotal::Integer,
    x_left::Real,
    x_right::Real,
)
    comm = Topologies.comm()

    ids = props.ids

    (; perm, copy) = commbuf
    recv_left, recv_right = commbuf.recvs

    LCS.@profile backend "particles/boundary/pack" begin
        send_left, send_right = pack_props_if_overflow!(xs, props, commbuf, backend, nvalid, x_left, x_right)
    end
    nsend_left = div(length(send_left), N_PROPERTIES)
    nsend_right = div(length(send_right), N_PROPERTIES)

    length(send_left) <= length(recv_left) || throw(
        ArgumentError("Left send buffer size $(length(send_left)) exceeds receive buffer size $(length(recv_left))")
    )
    length(send_right) <= length(recv_right) || throw(
        ArgumentError("Right send buffer size $(length(send_right)) exceeds receive buffer size $(length(recv_right))"),
    )

    LCS.@profile backend "particles/boundary/compact" begin
        compact_props!(props, commbuf, nvalid, backend)
    end

    LCS.@profile backend "particles/boundary/comm" begin
        recv_left, status = MPI.Sendrecv!(
            send_right, recv_left, comm, MPI.Status; dest=rank_right, source=rank_left, sendtag=axis, recvtag=axis
        )
        nrecv_left = div(MPI.Get_count(status, eltype(send_left)), N_PROPERTIES)
        recv_left = Utils.sliceshape(recv_left, nrecv_left, N_PROPERTIES)

        recv_right, status = MPI.Sendrecv!(
            send_left, recv_right, comm, MPI.Status; dest=rank_left, source=rank_right, sendtag=axis, recvtag=axis
        )
        nrecv_right = div(MPI.Get_count(status, eltype(send_right)), N_PROPERTIES)
        recv_right = Utils.sliceshape(recv_right, nrecv_right, N_PROPERTIES)
    end

    nsend = nsend_left + nsend_right
    nrecv = nrecv_left + nrecv_right
    nvalid - nsend + nrecv <= ntotal ||
        throw(ArgumentError("Particle count $(nvalid - nsend + nrecv) exceeds buffer capacity $ntotal"))

    nvalid = nvalid - nsend

    LCS.@profile backend "particles/boundary/unpack" begin
        indices_to_copy = (1:nrecv_left) .+ nvalid
        @views for (i, prop) in enumerate(Tuple(props))
            copyto!(prop[indices_to_copy], recv_left[:, i])
        end
        nvalid += nrecv_left

        indices_to_copy = (1:nrecv_right) .+ nvalid
        @views for (i, prop) in enumerate(Tuple(props))
            copyto!(prop[indices_to_copy], recv_right[:, i])
        end
        nvalid += nrecv_right
    end
    (; nvalid, nsend, nrecv)
end

function filter_indices!(
    indices::Property{<:Integer},
    mask::Property{<:Bool},
    scan::Property{<:Integer},
    nvalid::Integer,
    backend::KA.Backend,
)
    Utils.cumsum!(scan, mask, nvalid)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            if mask[i]
                indices[scan[i]] = i
            end
        end
    end

    nscan = if nvalid == 0
        0
    else
        Parallel.hostaccess(scan, nvalid)
    end
    view(indices, 1:nscan)
end

function pack_props_in_indices!(
    sendbuf::Property{<:Real}, indices::Property{<:Integer}, props::Properties, backend::KA.Backend
)
    ids, props... = Tuple(props)

    sendbuf_valid = Utils.sliceshape(sendbuf, length(indices), N_PROPERTIES)
    Parallel.foraxes(backend, (eachindex(indices),)) do i
        @inbounds begin
            sendbuf_valid[i, 1] = ids[indices[i]]
            for (j, prop) in enumerate(props)
                sendbuf_valid[i, j + 1] = prop[indices[i]]
            end
            ids[indices[i]] = Particles.INVALID
        end
    end

    # unwrap reshape
    parent(sendbuf_valid)
end

function pack_props_if_overflow!(
    xs::Property{<:Real},
    props::Properties,
    commbuf::CommBuffer,
    backend::KA.Backend,
    nvalid::Integer,
    x_left::Real,
    x_right::Real,
)
    send_left, send_right = commbuf.sends
    (; indices, scan, mask) = commbuf

    ids = props.ids

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x = xs[i]
            mask[i] = x < x_left && ids[i] != Particles.INVALID
        end
    end
    indices_left = filter_indices!(indices, mask, scan, nvalid, backend)
    send_left = pack_props_in_indices!(send_left, indices_left, props, backend)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x = xs[i]
            mask[i] = x > x_right && ids[i] != Particles.INVALID
        end
    end
    indices_right = filter_indices!(indices, mask, scan, nvalid, backend)
    send_right = pack_props_in_indices!(send_right, indices_right, props, backend)

    send_left, send_right
end

function compact_props!(props::Properties, commbuf::CommBuffer, nvalid::Integer, backend::KA.Backend)
    (; scan, mask, copy) = commbuf
    ids = props.ids

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds mask[i] = ids[i] != Particles.INVALID
    end

    Utils.cumsum!(scan, mask, nvalid)

    for prop in Tuple(props)
        c = reinterpret(eltype(prop), copy)
        copyto!(view(c, 1:nvalid), view(prop, 1:nvalid))
        KA.synchronize(backend)
        Parallel.foraxes(backend, (1:nvalid,)) do i
            @inbounds begin
                if mask[i]
                    prop[scan[i]] = c[i]
                end
            end
        end
    end
end
