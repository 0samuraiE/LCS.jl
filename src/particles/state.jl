"""
    ParticleState <: LCS.AbstractState

Particle-specific state information.
"""
mutable struct ParticleState <: LCS.AbstractState
    nvalid :: Int
end
@composite ParticleState

function ParticleState(::LCS.Config, ::Topologies.Topology)
    ParticleState(0)
end
