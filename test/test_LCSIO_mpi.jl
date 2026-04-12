module TestLCSIOMPI
using Test

using Accessors
using HDF5
using LCS
using MPI
using Offsets
using Random
using Topologies

import KernelAbstractions as KA

MPI.Initialized() || MPI.Init()

proc_dims = length(ARGS) == 3 ? parse.(Int, ARGS[1:3]) : [0, 0, 0]

backend = KA.CPU()
Topologies.device!(backend)

comm = Topologies.comm()
info = MPI.Info()

topo = Topologies.Topology(proc_dims)
rank = Topologies.linear_rank(topo)
cart_rank = Topologies.cart_rank(topo)
proc_dims = Topologies.proc_dims(topo)

config = LCSIO.readconfig(joinpath(@__DIR__, "../precompile.lcs-yaml"))

# low-level read/write (field)
let
    dims = (12, 12, 12)
    dims_l = div.(dims, proc_dims)
    A = reshape(1:prod(dims), dims...)
    indices = LCSIO.local_field_indices(dims_l, cart_rank)
    A_l = @views A[indices...]

    fname = MPI.bcast(tempname(), 0, comm)
    h5open(fname, "w", comm, info) do h
        LCSIO.write_field(h, "test", A_l, topo)
        A_read = similar(A_l)
        LCSIO.read_field(h, "test", A_read, topo)
        @test A_read == A_l
    end

    h5open(fname, "r", comm, info) do h
        @test h["test"][:, :, :] == A
    end
end

# low-level read/write (property)
let
    count_l = rand(1:12)
    count_g = MPI.Allreduce(count_l, MPI.SUM, comm)
    A = 1:count_g

    counts = MPI.Allgather(count_l, comm)
    indices = LCSIO.local_property_indices(counts, rank)
    A_l = @views A[indices]

    fname = MPI.bcast(tempname(), 0, comm)
    h5open(fname, "w", comm, info) do h
        LCSIO.write_property(h, "test", A_l, topo)
    end

    h5open(fname, "r", comm, info) do h
        @test h["test"][:] == A
        A_read = similar(A_l)
        LCSIO.read_property(h, "test", A_read, topo)
        @test A_read == A_l
    end
end

# flow fields
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    U, V, W = buf.flow.fields.Us
    P = buf.flow.fields.P

    rand!(U)
    rand!(V)
    rand!(W)
    rand!(P)

    U_orig = copy(U)
    V_orig = copy(V)
    W_orig = copy(W)
    P_orig = copy(P)

    file = joinpath(tmpdir, "test_flow.h5")
    LCSIO.save_flow_fields(file, buf, config, topo)

    fill!(U, 0.0)
    fill!(V, 0.0)
    fill!(W, 0.0)
    fill!(P, 0.0)

    LCSIO.load_flow_fields!(file, buf, state, config, topo)

    dims_l = LCS.dims_l(config.grid, topo)
    indices = CartesianIndices(map(n -> 1:n, dims_l))
    @offsetviews config.grid begin
        @test U[indices] ≈ U_orig[indices]
        @test V[indices] ≈ V_orig[indices]
        @test W[indices] ≈ W_orig[indices]
        @test P[indices] ≈ P_orig[indices]
    end
end

# particle counts
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    nvalid_orig = rand(1:16)
    state.particles[1].nvalid = nvalid_orig

    file = joinpath(tmpdir, "test_counts.h5")
    LCSIO.save_particle_counts(file, state, config, topo, 1)

    state.particles[1].nvalid = 0

    LCSIO.load_particle_counts!(file, buf, state, config, topo, 1)

    @test state.particles[1].nvalid == nvalid_orig
end

# particle position
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    nvalid = rand(10:16)
    state.particles[1].nvalid = nvalid

    xs, ys, zs = buf.particles[1].props.xss
    @views begin
        rand!(xs[1:nvalid])
        rand!(ys[1:nvalid])
        rand!(zs[1:nvalid])

        xs_orig = copy(xs[1:nvalid])
        ys_orig = copy(ys[1:nvalid])
        zs_orig = copy(zs[1:nvalid])
    end

    file = joinpath(tmpdir, "test_position.h5")
    LCSIO.save_particle_counts(file, state, config, topo, 1)
    LCSIO.save_particle_position(file, buf, state, config, topo, 1)

    fill!(xs, 0.0)
    fill!(ys, 0.0)
    fill!(zs, 0.0)

    LCSIO.load_particle_position!(file, buf, state, config, topo, 1)

    @test xs[1:nvalid] ≈ xs_orig
    @test ys[1:nvalid] ≈ ys_orig
    @test zs[1:nvalid] ≈ zs_orig
end

# particle velocity
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    nvalid = rand(10:16)
    state.particles[1].nvalid = nvalid

    us, vs, ws = buf.particles[1].props.uss
    @views begin
        rand!(us[1:nvalid])
        rand!(vs[1:nvalid])
        rand!(ws[1:nvalid])

        us_orig = copy(us[1:nvalid])
        vs_orig = copy(vs[1:nvalid])
        ws_orig = copy(ws[1:nvalid])
    end

    file = joinpath(tmpdir, "test_velocity.h5")
    LCSIO.save_particle_counts(file, state, config, topo, 1)
    LCSIO.save_particle_velocity(file, buf, state, config, topo, 1)

    fill!(us, 0.0)
    fill!(vs, 0.0)
    fill!(ws, 0.0)

    LCSIO.load_particle_velocity!(file, buf, state, config, topo, 1)

    @test us[1:nvalid] ≈ us_orig
    @test vs[1:nvalid] ≈ vs_orig
    @test ws[1:nvalid] ≈ ws_orig
end

# particle size
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    nvalid = rand(10:16)
    state.particles[1].nvalid = nvalid

    diams = buf.particles[1].props.diams

    @views begin
        rand!(diams[1:nvalid])

        diams_orig = copy(diams[1:nvalid])
    end

    file = joinpath(tmpdir, "test_size.h5")
    LCSIO.save_particle_counts(file, state, config, topo, 1)
    LCSIO.save_particle_size(file, buf, state, config, topo, 1)

    fill!(diams, 0.0)

    LCSIO.load_particle_size!(file, buf, state, config, topo, 1)

    @test diams[1:nvalid] ≈ diams_orig
end

# particle id
let config = config
    tmpdir = MPI.bcast(mktempdir(), 0, comm)
    config = @set config.outdir = tmpdir

    buf = LCS.Buffer(config, topo)
    state = LCS.State(config, topo)

    nvalid = rand(10:16)
    state.particles[1].nvalid = nvalid
    nvalid_g = Particles.nvalid_g(config.particles[1].population)

    ids = buf.particles[1].props.ids
    for i in 1:nvalid
        ids[i] = i + Topologies.linear_rank(topo) * nvalid_g
    end

    @views ids_orig = copy(ids[1:nvalid])

    file = joinpath(tmpdir, "test_id.h5")
    LCSIO.save_particle_counts(file, state, config, topo, 1)
    LCSIO.save_particle_id(file, buf, state, config, topo, 1)

    fill!(ids, 0)

    LCSIO.load_particle_id!(file, buf, state, config, topo, 1)

    @test @views ids[1:nvalid] == ids_orig
end

# upsample
let
    tmpdir_coarse = MPI.bcast(tempname(), 0, comm)
    dims_coarse = div.(config.grid.dims, 2)
    config_coarse = @set config.grid.dims = dims_coarse
    config_coarse = @set config_coarse.outdir = tmpdir_coarse

    buf_coarse = LCS.Buffer(config_coarse, topo)
    state_coarse = LCS.State(config_coarse, topo)

    U_c, V_c, W_c = buf_coarse.flow.fields.Us
    P_c = buf_coarse.flow.fields.P

    rand!(U_c)
    rand!(V_c)
    rand!(W_c)
    rand!(P_c)

    config_coarse = LCSIO.init(config_coarse, topo)
    file_coarse = joinpath(tmpdir_coarse, "restart.h5.200")
    LCSIO.saverestart(file_coarse, buf_coarse, state_coarse, config_coarse, topo)

    tmpdir_fine = MPI.bcast(mktempdir(), 0, comm)
    dims_fine = 2 .* dims_coarse
    config_fine = @set config.grid.dims = dims_fine
    config_fine = @set config_fine.outdir = tmpdir_fine

    buf_fine = LCS.Buffer(config_fine, topo)
    state_fine = LCS.State(config_fine, topo)

    LCSIO.load_coarse_flow_fields!(file_coarse, buf_fine, state_fine, config_fine, topo)

    U_fine, V_fine, W_fine = buf_fine.flow.fields.Us
    P_fine = buf_fine.flow.fields.P

    dims_l_fine = LCS.dims_l(config_fine.grid, topo)
    indices_fine = CartesianIndices(map(n -> 1:n, dims_l_fine))

    @offsetviews config_fine.grid begin
        @test !all(iszero, U_fine[indices_fine])
        @test !all(iszero, V_fine[indices_fine])
        @test !all(iszero, W_fine[indices_fine])
        @test !all(iszero, P_fine[indices_fine])
    end

    config_fine = @set config_fine.flow.init = Flows.Upsample(file_coarse)
    buf_fine = LCS.Buffer(config_fine, topo)
    state_fine = LCS.State(config_fine, topo)

    @test_nowarn Flows.init!(buf_fine, state_fine, config_fine, topo)
end
end
