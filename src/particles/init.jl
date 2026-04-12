"""
    init!(buf, state, config, topo, iprofile)

Initialize all particle properties for a specific particle set.
"""
function init!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer)
    (; id, position, velocity, size) = config.particles[iprofile].init
    init_state!(position, buf, state, config, topo, iprofile)
    init_id!(id, buf, state, config, topo, iprofile)
    init_position!(position, buf, state, config, topo, iprofile)
    init_velocity!(velocity, buf, state, config, topo, iprofile)
    init_size!(size, buf, state, config, topo, iprofile)
end

function init_state!(
    ::InitPosition, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    nvalid_l = Particles.nvalid_l(config.particles[iprofile].population, topo)
    state.particles[iprofile].nvalid = nvalid_l
end

function init_state!(
    position::ParticleRestart,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    LCSIO.load_particle_counts!(position.file, buf, state, config, topo, position.iprofile)
end

function init_id!(
    id::InitId, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, iprofile::Integer
)
    init_id!(id; kwargs_id(id, buf, state, config, topo, iprofile)...)
end

function id_offset(nvalid_l::Integer, topo::Topologies.Topology)
    id_offset(Topologies.allgather(nvalid_l, topo), Topologies.linear_rank(topo))
end
function id_offset(nvalid_l_all::Vector{<:Integer}, rank::Integer)
    sum(nvalid_l_all[1:rank]; init=0)
end

function kwargs_id(::GenerateId, buf, state, config, topo, iprofile)
    (; buf.particles[iprofile].props.ids, state.particles[iprofile].nvalid, config.backend, topo)
end

function init_id!(
    ::GenerateId; ids::Property{<:Integer}, nvalid::Integer, backend::KA.Backend, topo::Topologies.Topology
)
    offset = id_offset(nvalid, topo)
    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            ids[i] = i + offset
        end
    end
end

function init_id!(
    id::ParticleRestart,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    LCSIO.load_particle_id!(id.file, buf, state, config, topo, id.iprofile)
end

function init_position!(
    position::InitPosition,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    init_position!(position; kwargs_position(position, buf, state, config, topo, iprofile)...)
end

function kwargs_position(::RandomPosition, buf, state, config, topo, iprofile)
    (; buf.particles[iprofile].props.xss, state.particles[iprofile].nvalid, config.backend, topo)
end

function init_position!(
    ::RandomPosition; xss::Tuple3{Property{<:Real}}, nvalid::Integer, backend::KA.Backend, topo::Topologies.Topology
)
    xs_valid, ys_valid, zs_valid = view.(xss, Ref(1:nvalid))

    copyto!(xs_valid, rand(eltype(xs_valid), nvalid))
    copyto!(ys_valid, rand(eltype(ys_valid), nvalid))
    copyto!(zs_valid, rand(eltype(zs_valid), nvalid))

    ox, oy, oz = LCS.origins_l(topo)
    lx, ly, lz = LCS.lengths_l(topo)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x, y, z = xs_valid[i], ys_valid[i], zs_valid[i]

            xs_valid[i] = ox + x * lx
            ys_valid[i] = oy + y * ly
            zs_valid[i] = oz + z * lz
        end
    end
end

function init_position!(
    position::ParticleRestart,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    LCSIO.load_particle_position!(position.file, buf, state, config, topo, position.iprofile)
end

function init_velocity!(
    velocity::InitVelocity,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    init_velocity!(velocity; kwargs_velocity(velocity, buf, state, config, topo, iprofile)...)
end

function kwargs_velocity(::RestVelocity, buf, state, config, topo, iprofile)
    (; buf.particles[iprofile].props.uss, state.particles[iprofile].nvalid, config.backend)
end

function init_velocity!(::RestVelocity; uss::Tuple3{Property{<:Real}}, nvalid::Integer, backend::KA.Backend)
    us, vs, ws = uss

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            us[i] = 0.0
            vs[i] = 0.0
            ws[i] = 0.0
        end
    end
end

function init_velocity!(
    velocity::ParticleRestart,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    LCSIO.load_particle_velocity!(velocity.file, buf, state, config, topo, velocity.iprofile)
end

function init_size!(
    diameter::InitSize,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    init_size!(diameter; kwargs_size(diameter, buf, state, config, topo, iprofile)...)
end

function kwargs_size(::ConstSize, buf, state, config, topo, iprofile)
    (;
        buf.particles[iprofile].props.diams,
        state.particles[iprofile].nvalid,
        fparams=config.flow.params,
        config.backend,
    )
end

function init_size!(
    size::ConstSize; diams::Property{<:Real}, nvalid::Integer, fparams::Flows.FlowParams, backend::KA.Backend
)
    diam = Particles.diam(size, fparams)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            diams[i] = diam
        end
    end
end

function init_size!(
    size::ParticleRestart,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    iprofile::Integer,
)
    LCSIO.load_particle_size!(size.file, buf, state, config, topo, size.iprofile)
end
