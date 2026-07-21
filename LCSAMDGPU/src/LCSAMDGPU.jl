module LCSAMDGPU
using Accessors
using AMDGPU
using LCS
using MPI
using Topologies

MPI.Init()

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        MPI.Initialized() || MPI.Init()

        comm = Topologies.comm()
        rank = MPI.Comm_rank(comm)

        tmpdir = if rank == 0
            mktempdir()
        else
            ""
        end
        tmpdir = MPI.bcast(tmpdir, 0, comm)

        config = LCSIO.readconfig(joinpath(dirname(pathof(LCS)), "../precompile.lcs-yaml"))
        @reset config.backend = AMDGPU.ROCBackend()
        @reset config.flow.stat.integral_length = false
        Topologies.device!(config)

        ENV["LCS_LOG_LEVEL"] = "QUIET"
        cd(tmpdir) do
            LCS.simulate(config)
        end
        delete!(ENV, "LCS_LOG_LEVEL")

        if rank == 0 && isdir(tmpdir)
            rm(tmpdir; recursive=true)
        end

        if !MPI.Finalized()
            MPI.Finalize()
        end
    end
end
end
