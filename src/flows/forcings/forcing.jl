"""
    forcing!(buf, state, config, topo, stage)

Apply flow forcing based on the configured forcing method.
"""
function forcing!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage2)
    forcing!(config.flow.forcing, buf, state, config, topo, stage)
end

function forcing!(::NoForcing, ::LCS.Buffer, ::LCS.State, ::LCS.Config, ::Topologies.Topology, ::LCS.RKStage2)
    LCS.EMPTY_LOG
end
