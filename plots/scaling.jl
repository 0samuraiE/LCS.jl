using CairoMakie
using CSV

include("common.jl")
using .Common: MAKIE_FONTS_THEME, ticks, log10ticks, logbounds, _logrange

Makie.set_theme!(Attributes(; fonts=MAKIE_FONTS_THEME))

getitem(d, field) = getindex.(d, field)
filteritem(d, N, field) = getindex.(d[d.N .== N], field)

for device in ("cpu", "gpu")
    file = "data/$device-strong.csv"
    d = CSV.File(file)
    for col in (1, 2)
        for fontsize in (16, 24, 32)
            f = Figure(; size=(800 / col, 494), fontsize)
            ax = Axis(#
                f[1, 1];
                xlabel="Number of $(uppercase(device))s",
                ylabel="Wall Time (s)",
                yscale=log10,
                xscale=log10,
                xgridvisible=false,
                xticks=ticks(Int.(_logrange(1, 512, 2))),
                yticks=log10ticks(_logrange(0.01, 1000.0, 10)),
                ygridvisible=false,
                yminorticksvisible=true,
                yminorticks=IntervalsBetween(9),
                xticksize=16,
                yticksize=16,
                yminorticksize=8,
            )
            d_1500 = d[d.N .== 1500]
            scatterlines!(
                ax,
                filteritem(d, 1500, :devices),
                filteritem(d, 1500, :time);
                color=:blue,
                marker=:circle,
                markersize=fontsize,
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
                markersize=fontsize,
                label=L"N^3=3000^3,\:N_p=1500^3",
            )

            lines!(
                ax,
                filteritem(d, 3000, :devices),
                filteritem(d, 3000, :time)[1] ./
                (filteritem(d, 3000, :devices) ./ filteritem(d, 3000, :devices)[1]);
                color=:black,
                linestyle=:dash,
                label="Ideal",
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
    file = "data/$device-weak.csv"
    d = CSV.File(file)
    for col in (1, 2)
        for fontsize in (16, 24, 32)
            f = Figure(; size=(800 / col, 494), fontsize)
            ax = Axis(#
                f[1, 1];
                xlabel="Number of $(uppercase(device))s",
                ylabel="Wall Time (s)",
                yscale=log10,
                xscale=log10,
                xgridvisible=false,
                xticks=ticks((1:6) .^ 3),
                yticks=log10ticks(_logrange(0.01, 1000.0, 10)),
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
                markersize=fontsize,
                # label=L"N^3=750^3,\:N_p=325^3",
            )

            if device == "gpu"
                xs = getitem(d, :devices)
                ys = getitem(d, :time)

                lines!(
                    ax,
                    xs,
                    fill(ys[2], length(xs));
                    color=:black,
                    linestyle=:dash,
                    label="Ideal (pow-of-2)",
                )

                lines!(
                    ax,
                    xs[3:end],
                    fill(ys[3], length(xs) - 2);
                    color=:black,
                    linestyle=(:dashdot, :dense),
                    label="Ideal (non-pow-of-2)",
                )
            else
                hlines!(#
                    ax,
                    getitem(d, :time)[2];
                    color=:black,
                    linestyle=:dash,
                    label="Ideal",
                )
            end

            mmax = maximum(getitem(d, :time))
            mmin = minimum(getitem(d, :time))
            ymin, ymax = logbounds(mmin, mmax)
            ylims!(ax, 0.1 * ymin, ymax * 10)

            axislegend(ax)

            mkpath("figs")
            save("figs/$device-weak-$fontsize-$col.pdf", f)
        end
    end
end
