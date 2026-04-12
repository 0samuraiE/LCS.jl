module CUDAExt
using CUDA
using MPI
using Topologies

function Topologies.select_device!(::CUDABackend, comm::MPI.Comm, rank::Integer)
    lgrid = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank)
    lrank = MPI.Comm_rank(lgrid)
    if lrank < CUDA.ndevices()
        CUDA.device!(lrank)
    end
end
end
