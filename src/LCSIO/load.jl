function load_resume_state!(::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    rank = Topologies.linear_rank(topo)
    step, t = if rank == 0
        resume_file = find_latest_restart_file(config.outdir)
        step = extract_step(resume_file)

        stat_file = joinpath(config.outdir, FILE_STAT)

        t = h5open(stat_file, "r") do h
            read(h[string(step, "/", TAG_STATE, "/t")])
        end
        parse(Int, step), t
    else
        0, 0.0
    end

    step, t = MPI.Bcast((step, t), 0, Topologies.comm())

    state.step = step
    state.t = t
end

function load_flow_fields!(file::String, buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    U, V, W = buf.flow.fields.Us
    P = buf.flow.fields.P

    dims_l = LCS.dims_l(config.grid, topo)
    backend = config.backend

    _h5open(file, "r", topo) do h
        g = h[TAG_FLOW]
        domain = CartesianIndices(map(n -> 1:n, dims_l))
        @offsetviews config.grid begin
            @log backend "io/read/flow/U" read_field(g, TAG_FLOW_U, U[domain], topo)
            @log backend "io/read/flow/V" read_field(g, TAG_FLOW_V, V[domain], topo)
            @log backend "io/read/flow/W" read_field(g, TAG_FLOW_W, W[domain], topo)
            @log backend "io/read/flow/P" read_field(g, TAG_FLOW_P, P[domain], topo)
        end

        KA.synchronize(config.backend)
    end
end

function load_particle_properties!(
    file::String,
    props::Tuple{Vararg{Tuple{String,<:AbstractArray}}},
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    backend = config.backend
    proc_size = Topologies.proc_size(topo)

    _h5open(file, "r", topo) do h
        particle_group_name = string(TAG_PARTICLE, "/", iprofile)
        g = h[particle_group_name]
        length(g[TAG_PARTICLE_COUNTS]) == proc_size ||
            throw(ArgumentError("process count mismatch, expected $proc_size got $(length(g[TAG_PARTICLE_COUNTS]))"))

        nvalid = state.particles[iprofile].nvalid
        @views for (key, A) in props
            @log backend "io/read/particle/$iprofile/$key" read_property(g, key, A[1:nvalid], topo)
        end
    end
end

# Must be called before loading particle properties.
function load_particle_counts!(
    file::String, ::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    backend = config.backend
    rank = Topologies.linear_rank(topo)
    proc_size = Topologies.proc_size(topo)

    _h5open(file, "r", topo) do h
        gname = string(TAG_PARTICLE, "/", iprofile)
        g = h[gname]
        length(g[TAG_PARTICLE_COUNTS]) == proc_size ||
            throw(ArgumentError("process count mismatch, expected $proc_size got $(length(g[TAG_PARTICLE_COUNTS]))"))

        @log backend "io/read/particle/$iprofile/counts" begin
            state.particles[iprofile].nvalid = g[TAG_PARTICLE_COUNTS][rank + 1]
        end
    end
end

function load_particle_position!(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    xs, ys, zs = buf.particles[iprofile].props.xss
    props = ((TAG_PARTICLE_X, xs), (TAG_PARTICLE_Y, ys), (TAG_PARTICLE_Z, zs))
    load_particle_properties!(file, props, state, config, topo, iprofile)
end

function load_particle_velocity!(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    us, vs, ws = buf.particles[iprofile].props.uss
    props = ((TAG_PARTICLE_U, us), (TAG_PARTICLE_V, vs), (TAG_PARTICLE_W, ws))
    load_particle_properties!(file, props, state, config, topo, iprofile)
end

function load_particle_size!(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    diams = buf.particles[iprofile].props.diams
    props = ((TAG_PARTICLE_DIAM, diams),)
    load_particle_properties!(file, props, state, config, topo, iprofile)
end

function load_particle_id!(
    file::String, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    ids = buf.particles[iprofile].props.ids
    props = ((TAG_PARTICLE_ID, ids),)
    load_particle_properties!(file, props, state, config, topo, iprofile)
end
