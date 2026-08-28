using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, log10ticks, logbounds, lograngeby, ticks

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))
d = CSV.File("data/cpu-strong-intra.csv")

processes = collect(d.processes)
times = collect(d.time)
ideal = times[1] ./ processes

for fontsize in (18, 24)
    f = Figure(; size=(800, 494), fontsize)
    ax = Axis(
        f[1, 1];
        xlabel="Number of CPU Cores",
        ylabel="Wall Time (s)",
        xscale=log10,
        yscale=log10,
        xgridvisible=false,
        ygridvisible=false,
        xticks=ticks(processes),
        yticks=log10ticks(lograngeby(0.1, 100.0, 10)),
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
        xticksize=16,
        yticksize=16,
        yminorticksize=8,
    )

    scatterlines!(
        ax,
        processes,
        times;
        color=:blue,
        marker=:circle,
        markersize=fontsize / 1.2,
        label="Measured",
    )
    lines!(ax, processes, ideal; color=:black, linestyle=:dash, label="Ideal")

    ymin, ymax = logbounds(minimum(times), maximum(times))
    ylims!(ax, 0.1 * ymin, 10 * ymax)

    mkpath("figs")
    save("figs/cpu-strong-intra-$fontsize.pdf", f)
end
