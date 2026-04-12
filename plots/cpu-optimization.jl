using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, ticks, log10ticks, logbounds, _logrange

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/cpu-optimization.csv")

for fontsize in (16, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel="Number of Threads Per Process",
        ylabel="Speedup",
        xscale=log10,
        xgridvisible=false,
        xticks=ticks(Int.(_logrange(1, 128, 2))),
        ygridvisible=false,
        xticksize=16,
        yticksize=16,
    )

    xs = d.threads
    ys = d.time[1] ./ d.time
    scatterlines!(ax, xs, ys; color=:blue, marker=:circle, markersize=fontsize)

    ymin, ymax = extrema(ys)
    ylims!(ax, 0.9 * ymin, ymax * 1.1)

    mkpath("figs")
    save("figs/cpu-optimization-$fontsize.pdf", f)
end
