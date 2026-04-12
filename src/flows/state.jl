"""
    FlowState <: LCS.AbstractState

Flow-specific state information.
"""
mutable struct FlowState <: LCS.AbstractState end
@composite FlowState

function FlowState(::LCS.Config, topo::Topologies.Topology)
    FlowState()
end
