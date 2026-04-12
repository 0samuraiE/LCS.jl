function Topologies.device!(config::Config; multiprocessing::Bool=MPI.Initialized())
    Topologies.device!(config.backend; multiprocessing)
end

function Topologies.Topology(config::Config)
    Topologies.Topology(collect(config.topo.proc_dims); config.topo.reorder)
end

get_max_comm_arrays() = 4
function Topologies.HaloBuffer(config::Config, topo::Topologies.Topology)
    Topologies.HaloBuffer(
        config.backend, topo, LCS.FP, dims_l(config.grid, topo), config.grid.halo_size, get_max_comm_arrays()
    )
end

function Topologies.FullHalo(grid::Grid)
    Topologies.FullHalo(grid.halo_size)
end

function Topologies.FaceHalo(grid::Grid)
    Topologies.FaceHalo(grid.halo_size)
end

function Topologies.synchalo!(
    compute::Function,
    style::Topologies.HaloStyle,
    As::Tuple{Vararg{Field}},
    buf::Topologies.HaloBuffer,
    grid::Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    Topologies.synchalo!(compute, style, As, buf, LCS.dims_l(grid, topo), grid.halo_size, backend, topo)
end

function Topologies.synchalo!(
    style::Topologies.HaloStyle,
    As::Tuple{Vararg{Field}},
    buf::Topologies.HaloBuffer,
    grid::Grid,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    Topologies.synchalo!(style, As, buf, LCS.dims_l(grid, topo), grid.halo_size, backend, topo)
end
