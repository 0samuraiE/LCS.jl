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

        config = LCSIO.readconfig(joinpath(@__DIR__, "../precompile.lcs-yaml"))
        Topologies.device!(config)

        ENV[LCS.ENV_LCS_LOG_LEVEL] = LCS.LCS_LOG_QUIET
        cd(tmpdir) do
            simulate(config)
        end
        delete!(ENV, LCS.ENV_LCS_LOG_LEVEL)

        if rank == 0 && isdir(tmpdir)
            rm(tmpdir; recursive=true)
        end

        if !MPI.Finalized()
            MPI.Finalize()
        end
    end
end
