module TestIntegrationMPI
using Test

using Accessors
using HDF5
using LCS
using MPI
using Topologies

import KernelAbstractions as KA

MPI.Initialized() || MPI.Init()

proc_dims = length(ARGS) == 3 ? parse.(Int, ARGS[1:3]) : [0, 0, 0]

backend = KA.CPU()
Topologies.device!(backend)

comm = Topologies.comm()

topo = Topologies.Topology(proc_dims)
rank = Topologies.linear_rank(topo)

config = LCSIO.readconfig(joinpath(@__DIR__, "../precompile.lcs-yaml"))

function assert_particle_props_equal(actual_props, expected_props, nvalid)
    @test actual_props.ids[1:nvalid] == expected_props.ids[1:nvalid]
    @test actual_props.diams[1:nvalid] == expected_props.diams[1:nvalid]

    for (actual, expected) in zip(actual_props.xss, expected_props.xss)
        @test actual[1:nvalid] == expected[1:nvalid]
    end

    for (actual, expected) in zip(actual_props.uss, expected_props.uss)
        @test actual[1:nvalid] == expected[1:nvalid]
    end
end

tmpdir = if rank == 0
    mktempdir()
else
    ""
end
tmpdir = MPI.bcast(tmpdir, 0, comm)

cd(tmpdir) do
    let config = LCSIO.init(config, topo)
        buffer = LCS.Buffer(config, topo)
        state = LCS.State(config, topo)

        LCS.init!(buffer, state, config, topo)

        LCS.timestep!(buffer, state, config, topo)
        log = LCS.evolve!(buffer, state, config, topo)

        Topologies.isroot(topo) && LCS.printstate(state)
        Topologies.isroot(topo) && LCS.printlog(log)

        LCSIO.oninterval(state, config.simulate.interval_stat) do
            stat = LCS.stat(buffer, state, config, topo)
            LCSIO.savestat(config, state, topo; stat, log)
        end

        LCSIO.oninterval(state, config.simulate.interval_restart) do
            LCSIO.saverestart(buffer, state, config, topo)
        end

        # Generated files
        outdir = config.outdir
        restart_file = joinpath(outdir, "restart.h5.1")

        if Topologies.isroot(topo)
            @test isdir(outdir)

            generated_files = sort(readdir(outdir))
            @test "config.lcs-yaml" in generated_files
            @test "stat.h5" in generated_files

            restart_files = filter(f -> occursin(r"restart\.h5\.\d+", f), generated_files)
            @test sort(restart_files) == ["restart.h5.1"]

            h5open(joinpath(outdir, "stat.h5"), "r") do h
                @test haskey(h, "1")
                @test haskey(h["1"], "state")
                @test haskey(h["1"], "stat")
                @test haskey(h["1"], "log")

                saved_state = LCS.Utils.deepread(h, "1/state")
                @test saved_state["step"] == 1

                saved_stat = LCS.Utils.deepread(h, "1/stat")
                @test haskey(saved_stat, "state")
                @test haskey(saved_stat, "flow")
                @test haskey(saved_stat, "particles")

                @test saved_stat["state"]["step"] == 1
            end

            @test isfile(restart_file)

            h5open(restart_file, "r") do h
                @test haskey(h, "flow")
                @test haskey(h, "particles")

                flow_group = h["flow"]
                @test haskey(flow_group, "U")
                @test haskey(flow_group, "V")
                @test haskey(flow_group, "W")
                @test haskey(flow_group, "P")

                particles_group = h["particles"]
                @test haskey(particles_group, "1")

                particle_group = particles_group["1"]
                @test haskey(particle_group, "counts")
                @test haskey(particle_group, "id")
                @test haskey(particle_group, "x")
                @test haskey(particle_group, "y")
                @test haskey(particle_group, "z")
                @test haskey(particle_group, "u")
                @test haskey(particle_group, "v")
                @test haskey(particle_group, "w")
                @test haskey(particle_group, "diam")

                counts = read(particle_group["counts"])
                @test sum(counts) == config.particles[1].population.valid
            end
        end
        MPI.Barrier(comm)

        # Load and verify data consistency
        buffer2 = LCS.Buffer(config, topo)
        state2 = LCS.State(config, topo)

        LCSIO.load_flow_fields!(restart_file, buffer2, state2, config, topo)
        LCSIO.load_particle_counts!(restart_file, buffer2, state2, config, topo, 1)
        LCSIO.load_particle_id!(restart_file, buffer2, state2, config, topo, 1)
        LCSIO.load_particle_position!(restart_file, buffer2, state2, config, topo, 1)
        LCSIO.load_particle_velocity!(restart_file, buffer2, state2, config, topo, 1)
        LCSIO.load_particle_size!(restart_file, buffer2, state2, config, topo, 1)

        @test state2.particles[1].nvalid == state.particles[1].nvalid

        U1, V1, W1 = buffer.flow.fields.Us
        U2, V2, W2 = buffer2.flow.fields.Us
        P1 = buffer.flow.fields.P
        P2 = buffer2.flow.fields.P
        @test size(U1) == size(U2)
        @test size(V1) == size(V2)
        @test size(W1) == size(W2)

        nvalid = state.particles[1].nvalid
        assert_particle_props_equal(buffer2.particles[1].props, buffer.particles[1].props, nvalid)

        grid = config.grid
        @test LCS.unhalo(U1, grid) ≈ LCS.unhalo(U2, grid)
        @test LCS.unhalo(V1, grid) ≈ LCS.unhalo(V2, grid)
        @test LCS.unhalo(W1, grid) ≈ LCS.unhalo(W2, grid)
        @test LCS.unhalo(P1, grid) ≈ LCS.unhalo(P2, grid)

        # Resume from saved restart file
        ENV[LCS.ENV_LCS_RESUME] = "1"
        config_resume = LCSIO.readconfig(joinpath(config.outdir, "config.lcs-yaml"))

        @test LCSIO.isresumable(config_resume)

        Topologies.device!(config_resume)
        buffer_resume = LCS.Buffer(config_resume, topo)
        state_resume = LCS.State(config_resume, topo)

        config_resume = LCSIO.init(config_resume, topo)
        @test config_resume.resume == true

        LCS.init!(buffer_resume, state_resume, config_resume, topo)

        @test state_resume.step == state.step
        @test state_resume.t == state.t
        @test state_resume.particles[1].nvalid == state.particles[1].nvalid

        U_orig, V_orig, W_orig = buffer.flow.fields.Us
        U_resume, V_resume, W_resume = buffer_resume.flow.fields.Us
        P_orig = buffer.flow.fields.P
        P_resume = buffer_resume.flow.fields.P

        grid = config.grid
        @test LCS.unhalo(U_orig, grid) ≈ LCS.unhalo(U_resume, config_resume.grid)
        @test LCS.unhalo(V_orig, grid) ≈ LCS.unhalo(V_resume, config_resume.grid)
        @test LCS.unhalo(W_orig, grid) ≈ LCS.unhalo(W_resume, config_resume.grid)
        @test LCS.unhalo(P_orig, grid) ≈ LCS.unhalo(P_resume, config_resume.grid)

        nvalid_resume = state_resume.particles[1].nvalid
        assert_particle_props_equal(buffer_resume.particles[1].props, buffer.particles[1].props, nvalid_resume)

        # Continue from resumed state
        LCS.timestep!(buffer_resume, state_resume, config_resume, topo)
        log_resume = LCS.evolve!(buffer_resume, state_resume, config_resume, topo)

        @test state_resume.step == state.step + 1
        @test state_resume.t > state.t

        stat_resume = LCS.stat(buffer_resume, state_resume, config_resume, topo)
        LCSIO.savestat(config_resume, state_resume, topo; stat=stat_resume, log=log_resume)
        LCSIO.saverestart(buffer_resume, state_resume, config_resume, topo)

        if Topologies.isroot(topo)
            restart_file_2 = joinpath(outdir, "restart.h5.2")
            @test isfile(restart_file_2)

            h5open(joinpath(outdir, "stat.h5"), "r") do h
                @test haskey(h, "1")
                @test haskey(h, "2")
                @test haskey(h["2"], "state")
                @test haskey(h["2"], "stat")
                @test haskey(h["2"], "log")

                saved_state_2 = LCS.Utils.deepread(h, "2/state")
                @test saved_state_2["step"] == 2

                saved_stat_2 = LCS.Utils.deepread(h, "2/stat")
                @test saved_stat_2["state"]["step"] == 2
            end

            h5open(restart_file_2, "r") do h
                @test haskey(h, "flow")
                @test haskey(h, "particles")

                h5open(restart_file, "r") do h_first
                    @test Set(keys(h["flow"])) == Set(keys(h_first["flow"]))
                    @test Set(keys(h["particles"])) == Set(keys(h_first["particles"]))
                    @test Set(keys(h["particles/1"])) == Set(keys(h_first["particles/1"]))
                end
            end
        end
        MPI.Barrier(comm)

        delete!(ENV, LCS.ENV_LCS_RESUME)
    end
end

if Topologies.isroot(topo) && isdir(tmpdir)
    rm(tmpdir; recursive=true)
end
end
