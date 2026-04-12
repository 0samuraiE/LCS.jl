module TestFlowsMPI
using Test

using LCS: LCS, Utils
using MPI
using Offsets
using Random
using Topologies

import LCS.Flows as F
import KernelAbstractions as KA

MPI.Initialized() || MPI.Init()

proc_dims = length(ARGS) == 3 ? parse.(Int, ARGS[1:3]) : [0, 0, 0]

backend = KA.CPU()
Topologies.device!(backend)
topo = Topologies.Topology(proc_dims)

grid = LCS.Grid((24, 24, 24), 3)
forcing = F.EnergyPreserveRCF(1.0, 1, 1, 6)

fields = F.Fields(backend, grid, topo)
forcingbuf = F.RCFBuffer(forcing, backend, topo)

dims_l = LCS.dims_l(grid, topo)
filter_dims = F.filter_dims(forcing, grid)
coarse_dims_g = F.coarse_dims_g(forcing)

Us_coarse_g_ref = LCS.Utils.@ntuple(reshape(1:prod(coarse_dims_g), coarse_dims_g), LCS.N_DIMS)
Us_g = map(U_coarse_g -> repeat(U_coarse_g; inner=filter_dims), Us_coarse_g_ref)

domain = map(n -> (1:n), dims_l)
for (U, U_g) in zip(fields.Us, Us_g)
    @offsetviews grid U[domain...] .= $U_g[LCS.indices_l(grid, topo)...]
end

F.boxmean!(Topologies.processing(topo); fields.Us, forcingbuf, grid, forcing, backend, topo)
for (U_coarse_g, U_coarse_g_ref) in zip(forcingbuf.Us_coarse_g, Us_coarse_g_ref)
    @test U_coarse_g == U_coarse_g_ref
end
end
