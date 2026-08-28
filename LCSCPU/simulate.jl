using MPI

using Accessors
using LCS
using Topologies
using Topologies.Utils

function read_cpu_config(config_path, patch)
    mktemp() do cpu_config_path, io
        for line in eachline(config_path)
            println(io, replace(line, "cuda" => "cpu"))
        end
        close(io)
        return LCSIO.readconfig(cpu_config_path; patch)
    end
end

function main(args)
    isempty(args) && error("Usage: julia --project=LCSCPU LCSCPU/simulate.jl <config-file> [patch-expression]")

    config_path = first(args)
    patch_expression = length(args) >= 2 ? args[2] : ""
    patch = include_string(Main, "function (config)\n$patch_expression\nreturn config\nend")
    config = read_cpu_config(config_path, patch)

    @time LCS.simulate(config)
end

MPI.Initialized() || MPI.Init()
try
    main(ARGS)
finally
    MPI.Finalized() || MPI.Finalize()
end
