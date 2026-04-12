using Test

using MPI
using Topologies

import KernelAbstractions as KA

MPI.Initialized() || MPI.Init()

backend = KA.CPU()
Topologies.device!(backend)
topo = Topologies.Topology()

proc_dims = Topologies.proc_dims(topo)
linear_rank = Topologies.linear_rank(topo)

mocktopo = Topologies.Mock.mocktopology(; linear_rank, proc_dims)
@test topo == mocktopo

cart_rank = Topologies.cart_rank(topo)
dims_l = (2, 2, 2)
halo_size = 2

dims_g = proc_dims .* dims_l
dims_with_halo = dims_l .+ 2 .* halo_size
indices_domain = map(dims_l) do n
    (1:n) .+ halo_size
end

indices_l = map(cart_rank, dims_l) do ip, n
    (1:n) .+ ip * n
end
indices_with_halo = map(cart_rank, dims_l) do ip, nl
    ((1 - halo_size):(nl + halo_size)) .+ ip * nl
end
indices_with_halo_periodic = map(indices_with_halo, dims_g) do indice, ng
    mod1.(indice, ng)
end

A_g = collect(reshape(1:prod(dims_g), dims_g...))
As_l = (fill(0, dims_with_halo...), fill(0, dims_with_halo...))
for A_l in As_l
    @views A_l[indices_domain...] = A_g[indices_l...]
end

style = Topologies.FullHalo(halo_size)
buf = Topologies.HaloBuffer(backend, topo, Int, dims_l, halo_size, length(As_l))
Topologies.synchalo!(style, As_l, buf, dims_l, halo_size, backend, topo)

A_l_synched = A_g[indices_with_halo_periodic...]
for A_l in As_l
    @test A_l == A_l_synched
end

MPI.Finalize()
@test MPI.Finalized()
