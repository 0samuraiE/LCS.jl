module TestParticlesMPI
using Test

using LCS
using MPI
using Random
using Topologies

import LCS.Particles as P
import KernelAbstractions as KA

MPI.Initialized() || MPI.Init()

proc_dims = length(ARGS) == 3 ? parse.(Int, ARGS[1:3]) : [0, 0, 0]

backend = KA.CPU()
Topologies.device!(backend)
topo = Topologies.Topology(proc_dims)

population = P.Population(216, 10.0, 10.0)

props = P.Properties(backend, population, topo)
commbuf = P.CommBuffer(backend, population, topo)

ntotal_l = P.ntotal_l(population, topo)
nvalid_l = P.nvalid_l(population, topo)
nvalid_g = P.nvalid_g(population)

ids_g = 1:nvalid_g
ids_l = props.ids
@views ids_l[1:nvalid_l] .= ids_g[LCS.indices_l(nvalid_l, topo)]

# seed must be the same for all processes
rng = Random.MersenneTwister(1234)
xss_g = LCS.Utils.@ntuple(((rand(rng, LCS.FP, nvalid_g) .- 0.5) .* 1.1 .+ 0.5) .* LCS.DOMAIN_LENGTH, 3)
xss_g_periodic = map(xss_g) do xs_g
    LCS.Utils.wrap.(xs_g, LCS.DOMAIN_ORIGIN, LCS.DOMAIN_LENGTH)
end

xss_l = props.xss

@views for (xs_l, xs_g) in zip(xss_l, xss_g)
    xs_l[1:nvalid_l] .= xs_g[LCS.indices_l(nvalid_l, topo)]
end

nvalid_l = let nvalid = nvalid_l
    for i in 1:maximum(Topologies.proc_dims(topo))
        (; nvalid) = P._boundary!(; props.xss, props.xss2, props, commbuf, nvalid, population, backend, topo)

        @test Topologies.allsum(nvalid, topo) == nvalid_g

        @views for (xs_l, xs_g_periodic) in zip(xss_l, xss_g_periodic)
            @test xs_l[1:nvalid] == xs_g_periodic[ids_l[1:nvalid]]
        end
    end
    nvalid
end

@views for (xs_l, ox, ex) in zip(xss_l, LCS.origins_l(topo), LCS.ends_l(topo))
    @test all(ox .<= xs_l[1:nvalid_l] .< ex)
end
end
