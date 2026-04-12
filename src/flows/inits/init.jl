"""
    init!(buf, state, config, topo)

Initialize flow fields based on the configured initialization method.
"""
function init!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    init!(config.flow.init, buf, state, config, topo)
end

function init!(
    init::Union{Init,LCS.Restart}, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology
)
    init!(init; kwargs_init(init, buf, state, config, topo)...)
end

function init!(init::LCS.Restart, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    LCSIO.load_flow_fields!(init.file, buf, state, config, topo)
    LCS.EMPTY_LOG
end

function init!(init::Upsample, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    LCSIO.load_coarse_flow_fields!(init.file, buf, state, config, topo)
    LCS.EMPTY_LOG
end
