"""
    FlowStat <: LCS.AbstractStat

Flow statistics containing computed turbulence parameters.

# Fields
- `u`: Velocity statistics in x-direction
- `v`: Velocity statistics in y-direction
- `w`: Velocity statistics in z-direction
- `k`: Turbulent kinetic energy
- `ϵ`: Dissipation rate
- `τη`: Kolmogorov time scale
- `lμ`: Taylor microscale
- `lη`: Kolmogorov length scale
- `Reλ`: Taylor Reynolds number
- `urms`: Root mean square velocity
- `kη`: Kolmogorov wavenumber
- `skew`: Skewness of velocity gradient
- `flat`: Flatness of velocity gradient
- `L`: Integral length scale
- `Cϵ`: Dissipation coefficient (ϵL/u'³)
"""
@kwdef struct FlowStat <: LCS.AbstractStat
    #! format: off
    u    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    v    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    w    :: LCS.Summary{Float64} = LCS.Summary{Float64}()
    k    :: Float64 = 0.0
    ϵ    :: Float64 = 0.0
    τη   :: Float64 = 0.0
    lμ   :: Float64 = 0.0
    lη   :: Float64 = 0.0
    Reλ  :: Float64 = 0.0
    urms :: Float64 = 0.0
    kη   :: Float64 = 0.0
    skew :: Float64 = 0.0
    flat :: Float64 = 0.0
    L    :: Float64 = 0.0
    Cϵ   :: Float64 = 0.0
    #! format: on
end
@composite FlowStat
PPrint.PrintStyle(::Type{<:FlowStat}) = PPrint.Tree()

PPrint.rename(::Val{:urms}) = "u'"
PPrint.rename(::Val{:skew}) = "-S"
PPrint.rename(::Val{:flat}) = "F"

"""
    stat(buf, state, config, topo)

Compute flow statistics from current velocity field.
"""
function stat(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    stat(;
        buf.flow.fields.Us,
        buf.flow.scratches,
        statbuf=buf.flow.stat,
        halobuf=buf.halo,
        config.grid,
        config.flow.params,
        config.backend,
        topo,
    )
end
function stat(;
    Us::Tuple3{Field},
    scratches::Scratches,
    statbuf::StatBuffer,
    halobuf::Topologies.HaloBuffer,
    grid::LCS.Grid,
    params::FlowParams,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    U2, S2, D2, D3, D4 = scratches

    Topologies.synchalo!(Topologies.FullHalo(grid), Us, halobuf, grid, backend, topo)

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))
    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            local u = U[i, j, k]
            local v = V[i, j, k]
            local w = W[i, j, k]
            local u2 = u^2 + v^2 + w^2

            dudx = dxi * LCS.diff1_o4_at_half(U[(i - 2):(i + 1), j, k])
            dvdy = dyi * LCS.diff1_o4_at_half(V[i, (j - 2):(j + 1), k])
            dwdz = dzi * LCS.diff1_o4_at_half(W[i, j, (k - 2):(k + 1)])

            dudy = dyi * LCS.diff1_o4_at_half(U[i, (j - 2):(j + 1), k])
            dvdx = dxi * LCS.diff1_o4_at_half(V[(i - 1):(i + 2), j - 1, k])

            dvdz = dzi * LCS.diff1_o4_at_half(V[i, j - 1, (k - 1):(k + 2)])
            dwdy = dyi * LCS.diff1_o4_at_half(W[i, (j - 2):(j + 1), k])

            dwdx = dxi * LCS.diff1_o4_at_half(W[(i - 1):(i + 2), j, k - 1])
            dudz = dzi * LCS.diff1_o4_at_half(U[i, j, (k - 2):(k + 1)])

            s11 = 0.5 * (dudx + dudx)
            s22 = 0.5 * (dvdy + dvdy)
            s33 = 0.5 * (dwdz + dwdz)
            s12 = 0.5 * (dudy + dvdx)
            s13 = 0.5 * (dudz + dwdx)
            s23 = 0.5 * (dvdz + dwdy)

            local s2 = 2 * (s11^2 + s22^2 + s33^2 + 2 * (s12^2 + s13^2 + s23^2))

            S2[i, j, k] = s2
            U2[i, j, k] = u2
            D2[i, j, k] = dudx^2
            D3[i, j, k] = dudx^3
            D4[i, j, k] = dudx^4
        end
    end

    @offsetviews grid begin
        u = LCS.Summary(U[domain...], topo)
        v = LCS.Summary(V[domain...], topo)
        w = LCS.Summary(W[domain...], topo)

        u2 = Topologies.allmean(U2[domain...], topo)
        s2 = Topologies.allmean(S2[domain...], topo)

        d2 = Topologies.allmean(D2[domain...], topo)
        d3 = Topologies.allmean(D3[domain...], topo)
        d4 = Topologies.allmean(D4[domain...], topo)
    end

    urms = sqrt(u2 / 3)
    k = 0.5 * (u2 - u.mean^2 - v.mean^2 - w.mean^2)
    Re = Flows.Re(params)
    ϵ = s2 / Re
    τη = sqrt(1 / ϵ / Re)
    lμ = sqrt(15 * urms^2 / ϵ / Re)
    lη = (1 / ϵ / Re^3)^0.25
    Reλ = lμ * urms * Re

    n = minimum(LCS.dims_g(grid))
    kη = (n / 2) * lη

    skew = -d3 / d2^(3 / 2)
    flat = d4 / d2^2

    L = integral_length(U, statbuf.integral_length, grid, backend)
    Cϵ = ϵ * L / urms^3

    FlowStat(; u, v, w, k, ϵ, τη, lμ, lη, Reλ, urms, kη, skew, flat, L, Cϵ)
end

integral_length(::Field, ::Nothing, ::LCS.Grid, ::KA.Backend) = 0.0
function integral_length(U::Field, buf::IntegralLengthBuffer, grid::LCS.Grid, backend::KA.Backend)
    (; R, Ls, U_hat, plan) = buf

    R .= LCS.unhalo(U, grid)

    dx = LCS.spacings(grid)[1]
    n1, ns... = size(R)

    Utils.rfft!(U_hat, plan, R)
    U_hat .= abs2.(U_hat)
    Utils.irfft!(R, plan, U_hat)

    mid = div(n1, 2) + 1
    Parallel.foraxes(backend, map(n -> 1:n, ns)) do j, k
        @inbounds begin
            acc = 0.0
            for i in 1:mid
                r = R[i, j, k]
                if r < 0
                    break
                end
                acc += R[i, j, k]
            end
            Ls[j, k] = acc * dx / R[1, j, k]
        end
    end

    mean(Ls)
end

"""
    Q(buf, state, config, topo)

Compute the Q-criterion for vortex identification.
"""
function Q(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    Q(; buf.flow.fields.Us, buf.flow.scratches, halobuf=buf.halo, config.grid, config.backend, topo)
end
function Q(;
    Us::Tuple3{Field},
    scratches::Scratches,
    halobuf::Topologies.HaloBuffer,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    Uc, Vc, Wc, Out = scratches

    Topologies.synchalo!(Topologies.FullHalo(grid), Us, halobuf, grid, backend, topo)

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            u = LCS.interp_o4_at_half(U[(i - 2):(i + 1), j, k])
            v = LCS.interp_o4_at_half(V[i, (j - 2):(j + 1), k])
            w = LCS.interp_o4_at_half(W[i, j, (k - 2):(k + 1)])

            Uc[i, j, k] = u
            Vc[i, j, k] = v
            Wc[i, j, k] = w
        end
    end

    Topologies.synchalo!(Topologies.FullHalo(grid), (Uc, Vc, Wc), halobuf, grid, backend, topo)

    dxi, dyi, dzi = 1 ./ LCS.spacings(grid)

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            dudx = dxi * LCS.diff1_o4_at_node(Uc[(i - 2):(i + 2), j, k])
            dudy = dyi * LCS.diff1_o4_at_node(Uc[i, (j - 2):(j + 2), k])
            dudz = dzi * LCS.diff1_o4_at_node(Uc[i, j, (k - 2):(k + 2)])

            dvdx = dxi * LCS.diff1_o4_at_node(Vc[(i - 2):(i + 2), j, k])
            dvdy = dyi * LCS.diff1_o4_at_node(Vc[i, (j - 2):(j + 2), k])
            dvdz = dzi * LCS.diff1_o4_at_node(Vc[i, j, (k - 2):(k + 2)])

            dwdx = dxi * LCS.diff1_o4_at_node(Wc[(i - 2):(i + 2), j, k])
            dwdy = dyi * LCS.diff1_o4_at_node(Wc[i, (j - 2):(j + 2), k])
            dwdz = dzi * LCS.diff1_o4_at_node(Wc[i, j, (k - 2):(k + 2)])

            s2 = dudx^2 + dvdy^2 + dwdz^2 + 0.5 * ((dudy + dvdx)^2 + (dvdz + dwdy)^2 + (dwdx + dudz)^2)
            ω2 = 0.5 * ((dudy - dvdx)^2 + (dvdz - dwdy)^2 + (dwdx - dudz)^2)

            Out[i, j, k] = (ω2 - s2) / (ω2 + s2 + 1e-3)
        end
    end

    @offsetviews grid Out[domain...]
end

function speed(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    speed(; buf.flow.fields.Us, buf.flow.scratches, halobuf=buf.halo, config.grid, config.backend, topo)
end
function speed(;
    Us::Tuple3{Field},
    scratches::Scratches,
    halobuf::Topologies.HaloBuffer,
    grid::LCS.Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    Out, = scratches

    Topologies.synchalo!(Topologies.FullHalo(grid), Us, halobuf, grid, backend, topo)

    domain = map(n -> 1:n, LCS.dims_l(grid, topo))

    Parallel.foraxes(backend, domain) do i, j, k
        @inbounds @offsetviews grid begin
            u = LCS.interp_o4_at_half(U[(i - 2):(i + 1), j, k])
            v = LCS.interp_o4_at_half(V[i, (j - 2):(j + 1), k])
            w = LCS.interp_o4_at_half(W[i, j, (k - 2):(k + 1)])

            Out[i, j, k] = sqrt(u^2 + v^2 + w^2)
        end
    end

    @offsetviews grid Out[domain...]
end
