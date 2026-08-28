module Common

export KOLMOGOROV_COLORS, MAKIE_FONTS_THEME, log10ticks, logbounds, lograngeby, ticks

using Printf
using Makie

const MAKIE_FONTS_THEME = let
    MT = Makie.MathTeXEngine
    mt_fonts_dir = joinpath(
        dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern"
    )
    (;
        regular=joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
        bold=joinpath(mt_fonts_dir, "NewCM10-Bold.otf"),
    )
end

const KOLMOGOROV_COLORS = (:red, :blue, Makie.wong_colors()...)

ticks(x) = (x, string.(x))
function log10ticks(x)
    exps = ceil.(Int, log10.(x))
    labels = [L"10^{%$e}" for e in exps]
    return x, labels
end

function lograngeby(start, stop, base)
    length = ceil(Int, log(stop / start) / log(base)) + 1
    return logrange(start, stop; length)
end

function logbounds(mmin, mmax)
    ymin = 10.0^floor(log10(mmin))
    ymax = 10.0^ceil(log10(mmax))
    return ymin, ymax
end
end
