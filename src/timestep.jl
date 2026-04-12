"""
    timestep(buf, state, config, topo)

Compute adaptive time step satisfying CFL condition.
"""
function timestep(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    timestep(config.mode, buf, state, config, topo)
end

function timestep(::FlowMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    Flows.timestep(buf, state, config, topo)
end

function timestep(::FlowParticleMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    flow = Flows.timestep(buf, state, config, topo)
    particle = Particles.timestep(buf, state, config, topo)
    min(flow, particle)
end

"""
    timestep!(buf, state, config, topo)

Update `state.dt` with computed time step.
"""
function timestep!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    state.dt = timestep(buf, state, config, topo)
end
