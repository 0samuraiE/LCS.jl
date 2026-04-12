const SEPARATOR = "="^46

#! format: off
@kwdef struct Summary{T}
    mean :: Float64 = 0.0
    max  :: T = 0.0
    min  :: T = 0.0
end
#! format: on
@composite Summary

function PPrint.pprint(io::IO, s::Summary)
    PPrint.pprint(io, s.mean)
    print(io, "  ")
    PPrint.pprint(io, s.max)
    print(io, "  ")
    PPrint.pprint(io, s.min)
end

function Summary(A, topo::Topologies.Topology)
    Summary(Topologies.allmean(A, topo), Topologies.allmaximum(A, topo), Topologies.allminimum(A, topo))
end

abstract type AbstractStat end

struct StateStat <: AbstractStat
    step :: Int
    t    :: Float64
    dt   :: Float64
end
@composite StateStat
PPrint.PrintStyle(::Type{<:StateStat}) = PPrint.Tree()
StateStat(state::State) = StateStat(state.step, state.t, state.dt)

@concretize struct Stat <: AbstractStat
    state     :: StateStat
    flow      :: AbstractStat
    particles :: Tuple{Vararg{AbstractStat}}
end
@composite Stat
PPrint.PrintStyle(::Type{<:Stat}) = PPrint.Tree()

"""
    printstate(state)

Print simulation state (step, time, timestep).
"""
function printstate(state::LCS.State)
    t = @sprintf("step = %09d (t = %.2E, dt = %.2E)", state.step, state.t, state.dt)
    println(SEPARATOR)
    println(t)
    println(SEPARATOR)
end

"""
    printstat(stat; filter=identity)

Print simulation statistics with optional filtering.
"""
function printstat(stat; filter::Base.Callable=identity)
    label = "stat"
    header = string(SEPARATOR[1:3], " ", label, " ", SEPARATOR[(6 + length(label)):end])
    println(header)
    let io = IOContext(stdout, :realfmt => Printf.Format("%+.2E"))
        s = filter(stat)
        !isnothing(s) && PPrint.pprintln(io, s)
    end
    println(SEPARATOR)
end

"""
    printlog(log; filter=identity)

Print log data with optional filtering.
"""
function printlog(log; filter::Base.Callable=identity)
    label = "log"
    header = string(SEPARATOR[1:3], " ", label, " ", SEPARATOR[(6 + length(label)):end])
    println(header)
    let io = IOContext(stdout, :realfmt => Printf.Format("%+.2E"))
        l = filter(log)
        !isnothing(l) && PPrint.pprintln(io, l)
    end
    println(SEPARATOR)
end

"""
    stat(buf, state, config, topo)

Compute global statistics for current simulation state.
"""
function stat(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    stat(config.mode, buf, state, config, topo)
end

function stat(::FlowMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    flow = Flows.stat(buf, state, config, topo)
    Stat(StateStat(state), flow, ())
end

function stat(::FlowParticleMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    flow = Flows.stat(buf, state, config, topo)
    particles = ntuple(i -> Particles.stat(flow, buf, state, config, topo, i), length(buf.particles))
    Stat(StateStat(state), flow, particles)
end
