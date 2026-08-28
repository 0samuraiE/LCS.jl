using CairoMakie
using JLD2
using Printf

include("common.jl")
using .Common: MAKIE_FONTS_THEME, KOLMOGOROV_COLORS

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))

Ns = [128, 256, 512, 1024, 2048]
νs = [0.00699, 0.00278, 0.00110, 0.000450, 0.000179]
Reλs = [79.9, 129, 208, 333, 535]
lηs = [0.0328, 0.0160, 0.00783, 0.00416, 0.00211]
ϵs = [0.289, 0.330, 0.355, 0.311, 0.290]
for fontsize in (18, 24)
    f = Figure(; size=(750, 562), fontsize, figure_padding=(12, 24, 12, 12))
    ax = Axis(
        f[1, 1];
        xlabel=L"k\eta",
        ylabel=L"E(k) / (\epsilon \nu ^5)^{1/4}",
        xscale=log10,
        yscale=log10,
        xgridvisible=false,
        ygridvisible=false,
        xminorticks=IntervalsBetween(10),
        yminorticks=IntervalsBetween(10),
        xminorticksvisible=true,
        yminorticksvisible=true,
        xticksize=16,
        yticksize=16,
        xminorticksize=6,
        yminorticksize=6,
    )

    for (N, ν, Reλ, lη, ϵ, color) in zip(Ns, νs, Reλs, lηs, ϵs, KOLMOGOROV_COLORS)
        x = (1:(div(N, 2) + 1)) .* lη
        y = jldopen("data/energy-spectrum-$N.h5", "r")["data"]
        y = y / (ϵ * ν^5)^(1 / 4)

        lines!(ax, x, y; color, label=@sprintf("%.3g", Reλ))
    end
    x = lηs[2] .* (1:(Ns[2] ÷ 8))
    lines!(ax, x, 10^4.2 * x .^ (-5 / 3) / x[1]^(-5 / 3); color=:black)
    text!(ax, L"E(k) \sim k^{-5/3}"; position=(10^-1, 10^3))
    xlims!(ax, 1E-3, 1E1)
    ylims!(ax, 1E-4, 1E5)
    axislegend(ax, L"\mathrm{Re}_\lambda"; position=:rt, framevisible=false)

    mkpath("figs")
    save("figs/energy-spectrum-$fontsize.pdf", f)
end
