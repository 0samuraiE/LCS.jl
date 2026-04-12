Us_couple(::LCS.RKStage1, fields::Fields) = fields.Us2
Us_couple(::LCS.RKStage2, fields::Fields) = fields.Us

"""
    couple!(buf, state, config, topo, stage)

Apply pressure-velocity coupling to enforce incompressibility constraint.
"""
function couple!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage)
    couple!(config.flow.couple, buf, state, config, topo, stage)
end

function couple!(
    couple::Couple, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology, stage::LCS.RKStage
)
    couple!(couple; kwargs_couple(couple, buf, state, config, topo, stage)...)
end
