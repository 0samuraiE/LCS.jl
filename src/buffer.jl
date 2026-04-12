"""
    AbstractBuffer

Abstract base type for all buffer structs.
"""
abstract type AbstractBuffer end

"""
    Buffer <: AbstractBuffer

Top-level buffer struct containing all simulation buffers.
"""
@concretize struct Buffer <: AbstractBuffer
    halo      :: Topologies.HaloBuffer
    flow      :: AbstractBuffer
    particles :: Union{Tuple{Vararg{AbstractBuffer}},Nothing}
end

function Buffer(config::LCS.Config, topo::Topologies.Topology)
    Buffer(config.mode, config, topo)
end

function Buffer(::FlowMode, config::LCS.Config, topo::Topologies.Topology)
    Buffer(Topologies.HaloBuffer(config, topo), Flows.FlowBuffer(config, topo), nothing)
end

function Buffer(::FlowParticleMode, config::LCS.Config, topo::Topologies.Topology)
    particles = if !isnothing(config.particles)
        ntuple(i -> Particles.ParticleBuffer(config, topo, i), length(config.particles))
    else
        nothing
    end
    Buffer(Topologies.HaloBuffer(config, topo), Flows.FlowBuffer(config, topo), particles)
end
