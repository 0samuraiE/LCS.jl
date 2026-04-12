"""
    savestat(file, state, topo; stat, log=NamedTuple())
    savestat(config, state, topo; stat, log=NamedTuple())

Save simulation statistics and optional log data to HDF5 file.
"""
function savestat(
    config::LCS.Config, state::LCS.State, topo::Topologies.Topology; stat::LCS.Stat, log::NamedTuple=NamedTuple()
)
    file = joinpath(config.outdir, FILE_STAT)
    savestat(file, state, topo; stat, log)
end

function savestat(
    file::String, state::LCS.State, topo::Topologies.Topology; stat::LCS.Stat, log::NamedTuple=NamedTuple()
)
    Topologies.barrier(topo)

    if Topologies.isroot(topo)
        h5open(file, "cw") do h
            key = string(state.step)
            if haskey(h, key)
                delete_object(h, key)
            end
            g = create_group(h, key)
            Utils.deepwrite(g, TAG_STATE, PolySerde.normalize(state))
            Utils.deepwrite(g, TAG_STAT, PolySerde.normalize(stat))
            Utils.deepwrite(g, TAG_LOG, PolySerde.normalize(log))
        end
    end
    # Synchronize all ranks after root writes stat file
    Topologies.barrier(topo)
end

"""
    saverestart(buf, state, config, topo)

Save complete simulation state to HDF5 restart file.
"""
function saverestart(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    file = joinpath(config.outdir, Printf.format(FILE_RESTART_FMT, state.step))
    saverestart(file, buf, state, config, topo)
end

function saverestart(file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    save_flow_fields(file, buf, config, topo)
    if !isnothing(buf.particles)
        for i in eachindex(buf.particles)
            save_particle_counts(file, state, config, topo, i)
            save_particle_id(file, buf, state, config, topo, i)
            save_particle_position(file, buf, state, config, topo, i)
            save_particle_velocity(file, buf, state, config, topo, i)
            save_particle_size(file, buf, state, config, topo, i)
        end
    end
end

function save_flow_fields(file::String, buf::LCS.Buffer, config::LCS.Config, topo::Topologies.Topology)
    U, V, W = buf.flow.fields.Us
    P = buf.flow.fields.P

    dims_l = LCS.dims_l(config.grid, topo)
    backend = config.backend

    _h5open(file, "cw", topo) do h
        g = create_group(h, TAG_FLOW)

        domain = CartesianIndices(map(n -> 1:n, dims_l))
        # copying a non-contiguous GPU SubArray (view) cause internal GPU-side buffer allocation.
        A = Array(U)
        @offsetviews config.grid begin
            @log backend "io/write/flow/U" begin
                copyto!(A, U)
                write_field(g, TAG_FLOW_U, A[domain], topo)
            end
            @log backend "io/write/flow/V" begin
                copyto!(A, V)
                write_field(g, TAG_FLOW_V, A[domain], topo)
            end
            @log backend "io/write/flow/W" begin
                copyto!(A, W)
                write_field(g, TAG_FLOW_W, A[domain], topo)
            end
            @log backend "io/write/flow/P" begin
                copyto!(A, P)
                write_field(g, TAG_FLOW_P, A[domain], topo)
            end
        end
    end

    KA.synchronize(config.backend)
end

function save_particle_properties(
    file::String,
    props::Tuple{Vararg{Tuple{String,<:AbstractArray}}},
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    nvalid = state.particles[iprofile].nvalid
    backend = config.backend

    _h5open(file, "cw", topo) do h
        gname = string(TAG_PARTICLE, "/", iprofile)
        g = haskey(h, gname) ? h[gname] : create_group(h, gname)

        domain = 1:nvalid
        @views for (key, A) in props
            @log backend "io/write/particle/$iprofile/$key" write_property(g, key, A[domain], topo)
        end
    end

    KA.synchronize(config.backend)
end

function save_particle_counts(
    file::String, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    nvalid = state.particles[iprofile].nvalid
    backend = config.backend
    rank = Topologies.linear_rank(topo)
    proc_size = Topologies.proc_size(topo)

    _h5open(file, "cw", topo) do h
        gname = string(TAG_PARTICLE, "/", iprofile)
        g = haskey(h, gname) ? h[gname] : create_group(h, gname)
        @log backend "io/write/particle/$iprofile/counts" begin
            counts = _create_dataset(g, TAG_PARTICLE_COUNTS, Int64, (proc_size,), topo)
            # HDF5 uses 1-based indexing, but rank is 0-based, so add 1
            counts[rank + 1] = nvalid
        end
    end
end

function save_particle_position(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    xs, ys, zs = buf.particles[iprofile].props.xss
    props = ((TAG_PARTICLE_X, xs), (TAG_PARTICLE_Y, ys), (TAG_PARTICLE_Z, zs))
    save_particle_properties(file, props, state, config, topo, iprofile)
end

function save_particle_velocity(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    us, vs, ws = buf.particles[iprofile].props.uss
    props = ((TAG_PARTICLE_U, us), (TAG_PARTICLE_V, vs), (TAG_PARTICLE_W, ws))
    save_particle_properties(file, props, state, config, topo, iprofile)
end

function save_particle_size(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    diams = buf.particles[iprofile].props.diams
    props = ((TAG_PARTICLE_DIAM, diams),)
    save_particle_properties(file, props, state, config, topo, iprofile)
end

function save_particle_id(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    ids = buf.particles[iprofile].props.ids
    props = ((TAG_PARTICLE_ID, ids),)
    save_particle_properties(file, props, state, config, topo, iprofile)
end
