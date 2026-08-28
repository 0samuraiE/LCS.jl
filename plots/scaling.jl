using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, log10ticks, logbounds, lograngeby, ticks

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))

getitem(d, field) = getindex.(d, field)
filteritem(d, N, field) = getindex.(d[d.N .== N], field)

for device in ("cpu", "gpu")
    d = CSV.File("data/$device-strong.csv")
    for col in (1, 2)
        for fontsize in (18, 24)
            f = Figure(; size=(800 / col, 494), fontsize)
            ax = Axis(
                f[1, 1];
                xlabel=device == "cpu" ? "Number of CPU Sockets" : "Number of GPUs",
                ylabel="Wall Time (s)",
                yscale=log10,
                xscale=log10,
                xgridvisible=false,
                xticks=ticks(Int.(lograngeby(1, 512, 2))),
                yticks=log10ticks(lograngeby(0.01, 1000.0, 10)),
                ygridvisible=false,
                yminorticksvisible=true,
                yminorticks=IntervalsBetween(9),
                xticksize=16,
                yticksize=16,
                yminorticksize=8,
            )
            scatterlines!(
                ax,
                filteritem(d, 1500, :devices),
                filteritem(d, 1500, :time);
                color=:blue,
                marker=:circle,
                markersize=fontsize / 1.2,
                label=L"N^3=1500^3,\:N_p=750^3",
            )

            lines!(
                ax,
                filteritem(d, 1500, :devices),
                filteritem(d, 1500, :time)[1] ./
                (filteritem(d, 1500, :devices) ./ filteritem(d, 1500, :devices)[1]);
                color=:black,
                linestyle=:dash,
            )

            d_3000 = d[d.N .== 3000]
            scatterlines!(
                ax,
                getindex.(d_3000, :devices),
                getindex.(d_3000, :time);
                color=:red,
                marker=:xcross,
                markersize=fontsize / 1.2,
                label=L"N^3=3000^3,\:N_p=1500^3",
            )

            lines!(
                ax,
                filteritem(d, 3000, :devices),
                filteritem(d, 3000, :time)[1] ./
                (filteritem(d, 3000, :devices) ./ filteritem(d, 3000, :devices)[1]);
                color=:black,
                linestyle=:dash,
            )
            axislegend(ax; merge=true)

            mmax = max(
                maximum(filteritem(d, 1500, :time)), maximum(filteritem(d, 3000, :time))
            )
            mmin = min(
                minimum(filteritem(d, 1500, :time)), minimum(filteritem(d, 3000, :time))
            )
            ymin, ymax = logbounds(mmin, mmax)
            ylims!(ax, 0.1 * ymin, ymax * 10)

            mkpath("figs")
            save("figs/$device-strong-$fontsize-$col.pdf", f)
        end
    end
end

for device in ("cpu", "gpu")
    d = filter(row -> row.devices != 1, CSV.File("data/$device-weak.csv"))
    for col in (1, 2)
        for fontsize in (18, 24)
            f = Figure(; size=(800 / col, 494), fontsize)
            ax = Axis(
                f[1, 1];
                xlabel=device == "cpu" ? "Number of CPU Sockets" : "Number of GPUs",
                ylabel="Wall Time (s)",
                yscale=log10,
                xscale=log10,
                xgridvisible=false,
                xticks=ticks((2:6) .^ 3),
                yticks=log10ticks(lograngeby(0.01, 1000.0, 10)),
                ygridvisible=false,
                yminorticksvisible=true,
                yminorticks=IntervalsBetween(9),
                xticksize=16,
                yticksize=16,
                yminorticksize=8,
            )

            d = sort(d; by=e -> e.devices)

            scatterlines!(
                ax,
                getitem(d, :devices),
                getitem(d, :time);
                color=:blue,
                marker=:circle,
                markersize=fontsize / 1.2,
            )

            devices = getitem(d, :devices)
            baseline = getitem(d, :time)[only(findall(==(8), devices))]
            label = device == "gpu" ? "Ideal (8 GPUs)" : "Ideal (8 CPU Sockets)"
            hlines!(ax, baseline; color=:black, linestyle=:dash, label)

            mmax = maximum(getitem(d, :time))
            mmin = minimum(getitem(d, :time))
            ymin, ymax = logbounds(mmin, mmax)
            ylims!(ax, 0.1 * ymin, ymax * 10)

            mkpath("figs")
            save("figs/$device-weak-$fontsize-$col.pdf", f)
        end
    end
end
