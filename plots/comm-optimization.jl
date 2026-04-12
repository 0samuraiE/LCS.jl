using CairoMakie
using CSV

include("common.jl")
using .Common

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/comm-optimization.csv")

function _filter(d, method, field)
    return getindex.(d[d.method .== method], field)
end

for fontsize in (16, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel="Number of GPUs",
        ylabel="Speedup",
        xscale=log10,
        xgridvisible=false,
        xticks=ticks(Int.(_logrange(1, 128, 2))),
        ygridvisible=false,
        yminorticksvisible=true,
        xticksize=16,
        yticksize=16,
    )

    ymin, ymax = Inf, -Inf
    for (method, color, marker) in zip(
        ("overlap", "time-blocking", "time-blocking/overlap"),
        (:blue, :red, :green, nothing),
        (:circle, :xcross, :diamond, :utriangle),
    )
        xs = _filter(d, method, :devices)
        ys = _filter(d, "none", :time) ./ _filter(d, method, :time)

        scatterlines!(
            ax, xs, ys; color=color, marker=marker, markersize=fontsize, label=method
        )

        ymin = min(ymin, minimum(ys))
        ymax = max(ymax, maximum(ys))
    end

    ylims!(ax, 0.9 * ymin, ymax * 1.3)
    axislegend(ax)

    mkpath("figs")
    save("figs/comm-optimization-$fontsize.pdf", f)
end
