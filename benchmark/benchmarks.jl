using BenchmarkTools

using Accessors
using CUDA, cuDNN
using LCS
using Topologies

import KernelAbstractions as KA

const SUITE = BenchmarkGroup()

backends = Tuple{String,KA.Backend}[("cpu", KA.CPU())]
CUDA.has_cuda() && push!(backends, ("cuda", CUDA.CUDABackend()))

let
    for (name, backend) in backends
        let group = addgroup!(SUITE, name)
            config = LCSIO.readconfig(joinpath(@__DIR__, "benchmark.lcs-yaml"))
            config = @set config.backend = backend

            Topologies.device!(config)
            CUDA.device!(0)
            topo = Topologies.Topology(config)

            buffer = LCS.Buffer(config, topo)
            state = LCS.State(config, topo)

            let group = addgroup!(group, "flows")
                group["init"] = @benchmarkable begin
                    Flows.init!($buffer, $state, $config, $topo)
                    KA.synchronize($backend)
                end

                group["timestep"] = @benchmarkable begin
                    Flows.timestep!($buffer, $state, $config, $topo)
                    KA.synchronize($backend)
                end

                Flows.init!(buffer, state, config, topo)
                Flows.timestep!(buffer, state, config, topo)
                stage = LCS.RKStage2()

                group["rhs"] = @benchmarkable begin
                    Flows.rhs!($buffer, $state, $config, $topo, $stage)
                    KA.synchronize($backend)
                end
                group["update"] = @benchmarkable begin
                    Flows.update!($buffer, $state, $config, $topo, $stage)
                    KA.synchronize($backend)
                end
                group["couple"] = @benchmarkable begin
                    Flows.couple!($buffer, $state, $config, $topo, $stage)
                    KA.synchronize($backend)
                end
                group["forcing"] = @benchmarkable begin
                    Flows.forcing!($buffer, $state, $config, $topo, $stage)
                    KA.synchronize($backend)
                end

                group["stat"] = @benchmarkable begin
                    Flows.stat($buffer, $state, $config, $topo)
                    KA.synchronize($backend)
                end
            end

            if isnothing(buffer.particles)
                continue
            end
            iprofile = first(eachindex(buffer.particles))
            fstat = Flows.stat(buffer, state, config, topo)

            let group = addgroup!(group, "particles")
                group["init"] = @benchmarkable begin
                    Particles.init!($buffer, $state, $config, $topo, $iprofile)
                    KA.synchronize($backend)
                end

                group["timestep"] = @benchmarkable begin
                    Particles.timestep!($buffer, $state, $config, $topo)
                    KA.synchronize($backend)
                end

                Particles.init!(buffer, state, config, topo, iprofile)
                Particles.timestep!(buffer, state, config, topo)
                stage = LCS.RKStage2()

                group["makeindex"] = @benchmarkable begin
                    Particles.makeindex!($buffer, $state, $config, $topo, $stage, $iprofile)
                    KA.synchronize($backend)
                end

                group["motion"] = @benchmarkable begin
                    Particles.motion!($buffer, $state, $config, $topo, $stage, $iprofile)
                    KA.synchronize($backend)
                end

                group["update"] = @benchmarkable begin
                    Particles.update!($buffer, $state, $config, $topo, $stage, $iprofile)
                    KA.synchronize($backend)
                end

                group["stat"] = @benchmarkable begin
                    Particles.stat($fstat, $buffer, $state, $config, $topo, $iprofile)
                    KA.synchronize($backend)
                end
            end
        end
    end
end
