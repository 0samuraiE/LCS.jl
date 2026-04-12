using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, ticks, log10ticks, logbounds, _logrange

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/julia-vs-fortran.csv")

for fontsize in (16, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel="Number of Processes",
        ylabel="Wall Time (s)",
        yscale=log10,
        xscale=log10,
        xgridvisible=false,
        xticks=ticks([1, 2, 4, 8, 16, 32, 64]),
        yticks=log10ticks(logrange(0.1, 1000.0; length=5)),
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
        xticksize=16,
        yticksize=16,
        yminorticksize=8,
    )
    scatterlines!(
        ax, d.ps, d.julia; label="Julia", color=:blue, marker=:circle, markersize=fontsize
    )

    scatterlines!(
        ax,
        d.ps,
        d.fortran;
        label="Fortran",
        color=:red,
        marker=:xcross,
        markersize=fontsize,
    )
    axislegend(ax)
    ylims!(ax, 0.08, 1.25)

    mkpath("figs")
    save("figs/julia-vs-fortran-$fontsize.pdf", f)
end
