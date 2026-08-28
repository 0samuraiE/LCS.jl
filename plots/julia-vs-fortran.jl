using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, ticks

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/julia-vs-fortran.csv")
processes = collect(d.ps)

for fontsize in (16, 18, 20, 24, 32)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(#
        f[1, 1];
        xlabel="Number of CPU Cores",
        ylabel="Wall Time (s)",
        yscale=log10,
        xscale=log10,
        xgridvisible=false,
        xticks=ticks(processes),
        yticks=ticks([5, 10, 20, 30]),
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
        xticksize=16,
        yticksize=16,
        yminorticksize=8,
    )
    scatterlines!(
        ax,
        d.ps,
        d.julia;
        label="LCS.jl (Julia)",
        color=:blue,
        marker=:circle,
        markersize=fontsize / 1.2,
    )

    scatterlines!(
        ax,
        d.ps,
        d.fortran;
        label="LCS (Fortran)",
        color=:red,
        marker=:xcross,
        markersize=fontsize / 1.2,
    )
    axislegend(ax)
    ylims!(ax, 5, 35)

    mkpath("figs")
    save("figs/julia-vs-fortran-$fontsize.pdf", f)
end
