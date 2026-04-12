function rdf_assume_ci(;
    diam::AbstractFloat,
    xss::Tuple3{Property{<:Real}},
    props::Properties,
    cibuf::CellIndexBuffer,
    nvalid::Integer,
    grid::LCS.Grid,
    population::Population,
    cell::Cell,
    topo::Topologies.Topology,
    stat::Stat,
    backend::KA.Backend,
)
    xs, ys, zs = xss
    ids = props.ids
    (; hashes, starts, stops, hasher) = cibuf

    tol = stat.rdf.contact_tolerance
    nbins = stat.rdf.bin_size
    rmin = stat.rdf.min_radius
    rmax = stat.rdf.search_cell_size * minimum(LCS.spacings(grid))
    ncell = stat.rdf.search_cell_size
    edges = logrange(rmin, rmax, nbins + 1)

    npairs_d = KA.zeros(backend, Int, nbins)
    npairs_contact_d = KA.zeros(backend, Int, 1)

    @kernel function kernel!(::Val{NBINS}) where {NBINS}
        tid = @index(Global, Linear)
        lid = @index(Local, Linear)
        _npairs = @localmem Int NBINS
        _npairs_contact = @localmem Int 1

        gs = prod(@groupsize())
        k = lid
        while k <= NBINS
            _npairs[k] = 0
            k += gs
        end
        if lid == 1
            _npairs_contact[1] = 0
        end
        @synchronize()

        i1 = tid
        x1 = xs[i1], ys[i1], zs[i1]
        id1 = ids[i1]
        ci1 = Utils.decode(hasher, hashes[i1])

        if cell_in_domain(ci1, cell, grid, topo)
            for di in (-ncell):ncell, dj in (-ncell):ncell, dk in (-ncell):ncell
                ci2 = ci1 .+ (di, dj, dk)
                all(map((a, x) -> first(a) <= x <= last(a), hasher.axes, ci2)) || continue
                hash2 = Utils.encode(hasher, ci2)
                for i2 in starts[hash2]:stops[hash2]
                    id2 = ids[i2]
                    if id1 < id2
                        x2 = xs[i2], ys[i2], zs[i2]
                        d = Utils.norm(x1 .- x2)
                        k = Utils.searchsortedfirst(edges, eltype(edges)(d)) - 1
                        if 1 <= k <= NBINS
                            @atomic _npairs[k] += 1
                        end
                        if abs(d / diam - 1) < tol
                            @atomic _npairs_contact[1] += 1
                        end
                    end
                end
            end
        end
        @synchronize()

        gs = prod(@groupsize())
        k = lid
        while k <= NBINS
            @atomic npairs_d[k] += _npairs[k]
            k += gs
        end
        if lid == 1
            @atomic npairs_contact_d[1] += _npairs_contact[1]
        end
    end

    kernel!(backend)(Val(nbins); ndrange=nvalid)

    nvalid_g = Particles.nvalid_g(population)
    vol_domain_g = LCS.volume_g()
    rho0 = binomial(nvalid_g, 2) / vol_domain_g

    npairs_l = Parallel.hostaccess(npairs_d, 1:nbins)
    npairs = Topologies.allreduce(npairs_l, +, topo)
    gr = map(eachindex(npairs)) do k
        vol_shell = (4 / 3) * pi * (edges[k + 1]^3 - edges[k]^3)
        (npairs[k] / vol_shell) / rho0
    end

    npairs_contact_l = Parallel.hostaccess(npairs_contact_d, 1)
    npairs_contact = Topologies.allsum(npairs_contact_l, topo)
    vol_shell_contact = (4 / 3) * pi * diam^3 * 6 * tol
    gr_contact = (npairs_contact / vol_shell_contact) / rho0

    RDFStat(; gr, npairs, edges=collect(edges), gr_contact, npairs_contact)
end
