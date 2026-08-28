using CSV
using CairoMakie

include("common.jl")
using .Common: MAKIE_FONTS_THEME, KOLMOGOROV_COLORS, log10ticks, _logrange

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))

dat_ref = Dict(
    83 => Dict(0.2 => 3.94, 0.4 => 16.7, 0.6 => 26.2, 1.0 => 23.2, 2.0 => 7.59),
    133 => Dict(0.2 => 4.59, 0.4 => 16.6, 0.6 => 26.1, 1.0 => 24.8, 2.0 => 9.15),
    210 => Dict(0.2 => 4.67, 0.4 => 15.9, 0.6 => 25.7, 1.0 => 26.4, 2.0 => 10.8),
    346 => Dict(0.4 => 15.3, 0.6 => 24.3, 1.0 => 25.8, 2.0 => 13.1),
    537 => Dict(0.4 => 14.6, 0.6 => 23.2, 1.0 => 27.1, 2.0 => 14.7),
)

dat = Dict(
    80.2 => Dict(
        0.2 => 4.69,
        0.4 => 18.4,
        0.6 => 26.4,
        0.8 => 28.0,
        1.0 => 25.2,
        2.0 => 9.25,
        4.0 => 3.10,
        8.0 => 1.66,
    ),
    129.0 => Dict(
        0.2 => 5.11,
        0.4 => 16.5,
        0.6 => 26.1,
        0.8 => 28.6,
        1.0 => 26.5,
        2.0 => 11.6,
        4.0 => 3.95,
        8.0 => 2.07,
    ),
    209.0 => Dict(
        0.2 => 5.04,
        0.4 => 15.9,
        0.6 => 25.3,
        0.8 => 28.7,
        1.0 => 27.9,
        2.0 => 13.9,
        4.0 => 5.09,
        8.0 => 2.64,
    ),
    331.0 =>
        Dict(0.2 => 4.98, 0.4 => 15.2, 0.6 => 24.2, 0.8 => 28.0, 1.0 => 27.6, 2.0 => 15.0),
    534.0 => Dict(0.2 => 4.82, 0.4 => 14.4, 0.6 => 23.3, 1.0 => 27.9, 2.0 => 17.0),
)

params = CSV.File("data/parameterization.csv")[1]
model(St, Reλ) = model(St, Reλ, params)
function model(St, Reλ, params)
    A1 = params.A1
    A2 = params.A2
    ac = params.ac
    bc = params.bc
    ac2 = params.ac2
    bc2 = params.bc2
    aα = params.aα
    bα = params.bα

    y1 = A1 * St^2
    y2 = A2 * Reλ * St^(-2)
    St_a = (A2 / A1)^(1 / 4) * Reλ^(1 / 4)
    Ca = min(ac * Reλ^bc, ac2 * Reλ^bc2)
    α = log2(aα * Reλ^bα)
    za = 0.5 * (1 - tanh((log10(St) - log10(St_a)) / Ca))

    if St < St_a
        return 1 + y1 * za^α
    else
        return 1 + y2 * (1 - za)^α
    end
end

St_list = [0.4, 0.6, 1.0, 2.0]

for fontsize in (16, 18, 20, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel=L"\mathrm{Re}_\lambda",
        ylabel=L"g(r=R)",
        xscale=log10,
        yscale=identity,
        xgridvisible=false,
        ygridvisible=false,
        xticks=log10ticks(_logrange(0.01, 1000.0, 10)),
        xminorticks=IntervalsBetween(5),
        yminorticks=IntervalsBetween(5),
        xminorticksvisible=true,
        yminorticksvisible=true,
        xticksize=16,
        yticksize=16,
        xminorticksize=6,
        yminorticksize=6,
    )

    colors = KOLMOGOROV_COLORS
    for (St, color) in zip(St_list, colors)
        for Re in keys(dat_ref)
            scatter!(
                ax,
                Re,
                dat_ref[Re][St];
                marker=:circle,
                markersize=fontsize / 1.2,
                label="ref",
                color=:transparent,
                strokecolor=color,
                strokewidth=2,
            )
        end

        for Re in keys(dat)
            scatter!(
                ax,
                Re,
                dat[Re][St];
                color,
                marker=:x,
                markersize=fontsize / 1.2,
                label="present",
            )
        end

        x = logrange(10, 10000, 100)
        lines!(ax, x, model.(St, x); color)
    end

    axislegend(
        ax,
        [[
            MarkerElement(; color=c, marker=:circle, markersize=fontsize / 1.2) for
            c in colors[eachindex(St_list)]
        ],],
        [[L"%$St" for St in St_list]],
        [L"\mathrm{St}"];
        framevisible=true,
        position=:rb,
    )
    xlims!(ax, 50 / 1.1, 1000 * 1.1)
    ylims!(ax, 0, 35)
    f
    mkpath("figs")
    save("figs/gr-validation-$fontsize.pdf", f)
    fontsize == 24 && save("figs/gr-validation.pdf", f)
end
