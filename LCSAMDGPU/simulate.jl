using AMDGPU
using MPI

using Accessors
using LCS
using LCSAMDGPU

function main(args)
    isempty(args) && error("Usage: julia --project=LCSAMDGPU LCSAMDGPU/simulate.jl <config-file> [patch-expression]")

    config_path = first(args)
    patch_expression = length(args) >= 2 ? args[2] : ""
    patch = include_string(Main, "function (config)\n$patch_expression\nreturn config\nend")
    config = Base.invokelatest(LCSIO.readconfig, config_path; patch)

    LCS.simulate(config)
end

MPI.Initialized() || MPI.Init()
try
    main(ARGS)
finally
    MPI.Finalized() || MPI.Finalize()
end
