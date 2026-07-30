const WAVENUMBERS = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768)

function _fibonacci_direction(i::Integer, nsamples::Integer)
    phi = π * (3 - sqrt(5))
    y = 1 - (i / (nsamples - 1)) * 2
    radius = sqrt(1 - y * y)
    theta = phi * i
    cos(theta) * radius, y, sin(theta) * radius
end

function spectrum(;
    xss::Tuple3{Property{<:Real}},
    nvalid::Integer,
    population::Population,
    topo::Topologies.Topology,
    stat::Stat,
    backend::KA.Backend,
)
    nshells = Utils.searchsortedlast(WAVENUMBERS, stat.spectrum.max_wavenumber)
    nsamples = stat.spectrum.nsamples
    xs, ys, zs = xss
    nvalid_g = Particles.nvalid_g(population)
    cos_d = KA.zeros(backend, LCS.FP, nshells, nsamples)
    sin_d = KA.zeros(backend, LCS.FP, nshells, nsamples)

    @kernel function kernel!(::Val{NSHELLS}, ::Val{NSAMPLES}) where {NSHELLS,NSAMPLES}
        tid = @index(Global, Linear)
        lid = @index(Local, Linear)
        cos_local = @localmem LCS.FP (NSHELLS, NSAMPLES)
        sin_local = @localmem LCS.FP (NSHELLS, NSAMPLES)

        if lid == 1
            for b in 1:NSHELLS, i in 1:NSAMPLES
                cos_local[b, i] = 0
                sin_local[b, i] = 0
            end
        end
        @synchronize()

        x = xs[tid], ys[tid], zs[tid]
        for b in 1:NSHELLS, i in 1:NSAMPLES
            kk = WAVENUMBERS[b]
            k = kk .* _fibonacci_direction(i - 1, NSAMPLES)
            k_dot_x = sum(k .* x)
            @atomic cos_local[b, i] += cos(k_dot_x)
            @atomic sin_local[b, i] += sin(k_dot_x)
        end
        @synchronize()

        if lid == 1
            for b in 1:NSHELLS, i in 1:NSAMPLES
                @atomic cos_d[b, i] += cos_local[b, i]
                @atomic sin_d[b, i] += sin_local[b, i]
            end
        end
    end

    kernel!(backend)(Val(nshells), Val(nsamples); ndrange=nvalid)

    cos_l = Parallel.hostaccess(cos_d, 1:nshells, 1:nsamples)
    sin_l = Parallel.hostaccess(sin_d, 1:nshells, 1:nsamples)
    cos_g = Topologies.allreduce(cos_l, +, topo)
    sin_g = Topologies.allreduce(sin_l, +, topo)

    Enp = zeros(nshells)
    for b in 1:nshells
        phi = 0.0
        for i in 1:nsamples
            c = cos_g[b, i]
            s = sin_g[b, i]
            phi += (c^2 + s^2) / nvalid_g^2 - 1 / nvalid_g
        end
        kk = WAVENUMBERS[b]
        Enp[b] = phi / nsamples * 4 * π * kk^2
    end

    SpectrumStat(; wavenumbers=collect(WAVENUMBERS[1:nshells]), Enp)
end
