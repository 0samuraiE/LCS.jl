using Test

using Accessors
using LCS
using MPI
using Topologies

MPI.Initialized() || MPI.Init()

@testset "LCSIO" begin
    @testset "local indices" begin
        counts = [1, 2, 3, 4]
        @test LCSIO.local_property_indices(counts, 0) == 1:1
        @test LCSIO.local_property_indices(counts, 1) == 2:3
        @test LCSIO.local_property_indices(counts, 2) == 4:6
        @test LCSIO.local_property_indices(counts, 3) == 7:10

        @test LCSIO.local_field_indices((4, 4, 6), (0, 0, 0)) == (1:4, 1:4, 1:6)
        @test LCSIO.local_field_indices((4, 4, 6), (1, 0, 0)) == (5:8, 1:4, 1:6)
    end

    @testset "readconfig sync flag" begin
        cfg = LCSIO.readconfig(joinpath(@__DIR__, "../precompile.lcs-yaml"); sync=false)
        @test cfg isa LCS.Config

        tmpdir = mktempdir()
        source = joinpath(@__DIR__, "../precompile.lcs-yaml")
        dest = joinpath(tmpdir, "resume.lcs-yaml")
        open(dest, "w") do io
            write(io, read(source, String))
            write(io, "\nresume: true\n")
        end
        err = @test_throws ArgumentError LCSIO.readconfig(dest; sync=false)
        @test occursin("resume mode is controlled by LCS_RESUME environment variable", err.value.msg)
    end

    @testset "readconfig patch parameter" begin
        source = joinpath(@__DIR__, "../precompile.lcs-yaml")

        @testset "patch modifies config" begin
            cfg_default = LCSIO.readconfig(source; sync=false)
            original_cfl = cfg_default.timestep.cfl

            patch = config -> @set config.timestep.cfl = original_cfl * 2
            cfg_patched = LCSIO.readconfig(source; sync=false, patch)

            @test cfg_patched.timestep.cfl == original_cfl * 2
            @test cfg_patched.timestep.cfl != cfg_default.timestep.cfl
        end

        @testset "patch can modify outdir" begin
            tmpdir = joinpath(tempdir(), "test-patch-outdir-$(rand(UInt32))")
            patch = config -> @set config.outdir = tmpdir
            cfg = LCSIO.readconfig(source; sync=false, patch)

            @test cfg.outdir == tmpdir
        end

        @testset "patch applied before init" begin
            patch = config -> @set config.grid.dims = (32, 32, 32)
            cfg = LCSIO.readconfig(source; sync=false, patch)

            @test cfg.grid.dims == (32, 32, 32)
        end
    end

    @testset "isresumable detection" begin
        cfg = LCSIO.readconfig(joinpath(@__DIR__, "../precompile.lcs-yaml"); sync=false)
        tmpdir = mktempdir()
        cfg = @set cfg.outdir = tmpdir

        @test !LCSIO.isresumable(cfg)

        touch(joinpath(tmpdir, LCSIO.FILE_CONFIG))
        @test !LCSIO.isresumable(cfg)

        touch(joinpath(tmpdir, "restart.h5.1"))
        @test LCSIO.isresumable(cfg)
    end

    topo = Topologies.Topology([0, 0, 0])

    @testset "init" begin
        source = joinpath(@__DIR__, "../precompile.lcs-yaml")

        @testset "outdir exists without LCS_RESUME" begin
            tmpdir = mktempdir()
            cfg = LCSIO.readconfig(source; sync=false, patch=c -> @set c.outdir = tmpdir)

            err = @test_throws ArgumentError LCSIO.init(cfg, topo)
            @test occursin("already exists", err.value.msg)
            @test occursin("LCS_RESUME=1", err.value.msg)
        end

        @testset "outdir exists with LCS_RESUME but not resumable" begin
            tmpdir = mktempdir()
            cfg = LCSIO.readconfig(source; sync=false, patch=c -> @set c.outdir = tmpdir)

            ENV[LCS.ENV_LCS_RESUME] = "1"
            err = @test_throws ArgumentError LCSIO.init(cfg, topo)
            @test occursin("not resumable", err.value.msg)
            delete!(ENV, LCS.ENV_LCS_RESUME)
        end

        @testset "outdir exists with LCS_RESUME and resumable" begin
            tmpdir = mktempdir()
            touch(joinpath(tmpdir, LCSIO.FILE_CONFIG))
            touch(joinpath(tmpdir, "restart.h5.1"))
            cfg = LCSIO.readconfig(source; sync=false, patch=c -> @set c.outdir = tmpdir)

            ENV[LCS.ENV_LCS_RESUME] = "1"
            resolved = LCSIO.init(cfg, topo)
            @test resolved.resume == true
            delete!(ENV, LCS.ENV_LCS_RESUME)
        end

        @testset "outdir does not exist with LCS_RESUME" begin
            tmpdir = mktempdir()
            rm(tmpdir)
            cfg = LCSIO.readconfig(source; sync=false, patch=c -> @set c.outdir = tmpdir)

            ENV[LCS.ENV_LCS_RESUME] = "1"
            resolved = LCSIO.init(cfg, topo)
            @test resolved.resume == false
            @test resolved.outdir == tmpdir
            delete!(ENV, LCS.ENV_LCS_RESUME)
        end
    end
end
