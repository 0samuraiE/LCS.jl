abstract type AbstractState end

"""
    State

Simulation state container.

# Fields
- `step`: Current step number
- `t`: Current time
- `dt`: Current time step size
- `flow`: Flow state (nothing if not used)
- `particle`: Particle state (nothing if not used)
"""
mutable struct State{FlowT<:Union{AbstractState,Nothing},ParticlesT<:Union{Tuple{Vararg{AbstractState}},Nothing}} <:
               AbstractState
    step            :: Int
    t               :: Float64
    dt              :: Float64

    const flow      :: FlowT
    const particles :: ParticlesT
end
@composite State
PPrint.PrintStyle(::Type{<:State}) = PPrint.Tree()

function State(config::Config, topo::Topologies.Topology)
    flow = Flows.FlowState(config, topo)
    particles = if !isnothing(config.particles)
        ntuple(_ -> Particles.ParticleState(config, topo), length(config.particles))
    else
        nothing
    end
    State(0, 0.0, 0.0, flow, particles)
end

dtrk(::LCS.RKStage1, state::State) = 0.5 * state.dt
dtrk(::LCS.RKStage2, state::State) = state.dt
