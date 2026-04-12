"""
    init!(buf, state, config, topo)

Initialize simulation state and fields.
"""
function init!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    init_state!(buf, state, config, topo)
    init!(config.mode, buf, state, config, topo)
end

function init_state!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    if config.resume
        LCSIO.load_resume_state!(buf, state, config, topo)
    end
end

function init!(::FlowMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    Flows.init!(buf, state, config, topo)
end

function init!(::FlowParticleMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    Flows.init!(buf, state, config, topo)
    for i in eachindex(buf.particles)
        Particles.init!(buf, state, config, topo, i)
    end
end
