"""
    timestep(buf, state, config, topo)

Compute adaptive time step size for particle simulation.
"""
function timestep(buf::LCS.Buffer, ::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    if isnothing(buf.particles)
        return config.timestep.dtmax
    end
    minimum(eachindex(buf.particles)) do i
        timestep(buf.particles[i].props.diams, config.timestep, config.flow.params, config.particles[i].params, topo)
    end
end

function timestep(
    diams::Property{<:Real},
    ts::LCS.TimeStep,
    fprops::Flows.FlowParams,
    pprops::Particles.ParticleParams,
    topo::Topologies.Topology,
)
    diam = Topologies.allmaximum(diams, topo)
    taup = Particles.taup(pprops, fprops, diam)
    min(0.5 * taup, ts.dtmax)
end

"""
    timestep!(buf, state, config, topo)

Compute adaptive time step for particle simulation.
"""
function timestep!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    timestep(buf, state, config, topo)
end
