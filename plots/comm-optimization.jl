using CairoMakie
using CSV

include("common.jl")
using .Common

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/comm-optimization.csv")

function _filter(d, method, field)
    return getindex.(d[d.method .== method], field)
end

for fontsize in (16, 18, 20, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel="Number of GPUs",
        ylabel="Wall-Time Reduction (%)",
        xscale=log10,
        xgridvisible=false,
        xticks=ticks(Int.(_logrange(1, 128, 2))),
        ygridvisible=false,
        yminorticksvisible=true,
        xticksize=16,
        yticksize=16,
    )

    for (method, color, marker) in zip(
        ("overlap", "time-blocking", "time-blocking/overlap"),
        (:blue, :red, :green, nothing),
        (:circle, :xcross, :diamond, :utriangle),
    )
        xs = _filter(d, method, :devices)
        baseline = _filter(d, "none", :time)
        ys = 100 .* (baseline .- _filter(d, method, :time)) ./ baseline

        scatterlines!(
            ax, xs, ys; color=color, marker=marker, markersize=fontsize / 1.2, label=method
        )
    end

    ylims!(ax, 0, 15)
    axislegend(ax; position=:lt)

    mkpath("figs")
    save("figs/comm-optimization-$fontsize.pdf", f)
end
