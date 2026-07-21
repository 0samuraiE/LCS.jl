using CUDA, cuDNN
using MPI
MPI.Init()

using Accessors
using LCS
using LCSCUDA

length(ARGS) >= 1 || throw(ErrorException("Usage: julia simulate.jl <config file> [expression]"))
file = ARGS[1]
expr = length(ARGS) >= 2 ? ARGS[2] : ""

patch = include_string(Main, string("function (config)\n", expr, "\nreturn config\nend"))
config = Base.invokelatest(LCSIO.readconfig, file; patch)

LCS.simulate(config)

MPI.Finalize()
