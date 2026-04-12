using MPI
MPI.Init()

using Accessors
using LCS
using Topologies
using Topologies.Utils

length(ARGS) >= 1 || throw(ErrorException("Usage: julia simulate.jl <config file> [expression]"))
file = ARGS[1]
expr = length(ARGS) >= 2 ? ARGS[2] : ""

file_cpu, f = mktemp()
for line in eachline(file)
    if occursin("cuda", line)
        write(f, replace(line, "cuda" => "cpu"), "\n")
    else
        write(f, line, "\n")
    end
end
close(f)

patch = include_string(Main, string("function (config)\n", expr, "\nreturn config\nend"))

config = LCSIO.readconfig(file_cpu; patch)

@time LCS.simulate(config)

MPI.Finalize()

#=
launch commands examples:

mpirun -n 16 --map-by ppr:8:node \
    --bind-to core \
    --report-bindings
    -x LD_LIBRARY_PATH \
    julia --threads=1 --project=LCSCPU \
    LCSCPU/simulate.jl config/N128.lcs-yaml '@reset config.outdir="auto""'

mpirun -n 4 \
    --map-by ppr:1:package:pe=96 \
    --bind-to core \
    --report-bindings \
    -x LD_LIBRARY_PATH \
    julia --project -O3 --threads=96 \
    share/bench.jl share/config/N1500-bench.lcs-yaml \
    '@reset config.outdir *= "_4_96""'
=#
