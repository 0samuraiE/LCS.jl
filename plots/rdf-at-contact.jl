using CSV
using CairoMakie

include("common.jl")
using .Common: KOLMOGOROV_COLORS, MAKIE_FONTS_THEME, log10ticks, lograngeby

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))

rdf = CSV.File("data/rdf-at-contact.csv")
references = filter(row -> row.source == "reference", rdf)
measurements = filter(row -> row.source == "present", rdf)

params = CSV.File("data/rdf-model-parameters.csv")[1]
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

for fontsize in (18, 24)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(
        f[1, 1];
        xlabel=L"\mathrm{Re}_\lambda",
        ylabel=L"g(r=R)",
        xscale=log10,
        yscale=identity,
        xgridvisible=false,
        ygridvisible=false,
        xticks=log10ticks(lograngeby(0.01, 1000.0, 10)),
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
        for row in filter(row -> row.St == St, references)
            scatter!(
                ax,
                row.Re_lambda,
                row.g,
                marker=:circle,
                markersize=fontsize / 1.2,
                label="ref",
                color=:transparent,
                strokecolor=color,
                strokewidth=2,
            )
        end

        for row in filter(row -> row.St == St, measurements)
            scatter!(
                ax,
                row.Re_lambda,
                row.g;
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
    mkpath("figs")
    save("figs/rdf-at-contact-$fontsize.pdf", f)
end
