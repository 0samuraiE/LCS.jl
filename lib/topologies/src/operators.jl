"""
    allsum(f, A, topo::Topology)

Global sum reduction across all processes.
"""
allsum(f, A, topo::Topology) = _greduce(processing(topo), +, sum, f, A)
allsum(A, topo::Topology) = allsum(identity, A, topo)

"""
    allminimum(f, A, topo::Topology)

Global minimum reduction across all processes.
"""
allminimum(f, A, topo::Topology) = _greduce(processing(topo), min, minimum, f, A)
allminimum(A, topo::Topology) = allminimum(identity, A, topo)

"""
    allmaximum(f, A, topo::Topology)

Global maximum reduction across all processes.
"""
allmaximum(f, A, topo::Topology) = _greduce(processing(topo), max, maximum, f, A)
allmaximum(A, topo::Topology) = allmaximum(identity, A, topo)

"""
    allmean(f, A, topo::Topology)

Global mean reduction across all processes.
"""
allmean(f, A, topo::Topology) = _gmean_impl(processing(topo), f, A)
allmean(A, topo::Topology) = allmean(identity, A, topo)

@inline function _greduce(::SingleProcessing, op, reduce, f, A)
    reduce(f, A)
end

@inline function _greduce(::MultiProcessing, op, reduce, f, A)
    comm = Topologies.comm()
    MPI.Barrier(comm)
    lval = reduce(f, A)
    MPI.Allreduce(lval, op, comm)
end

@inline function _gmean_impl(::SingleProcessing, f, A)
    sum(f, A) / length(A)
end

@inline function _gmean_impl(::MultiProcessing, f, A)
    comm = Topologies.comm()
    MPI.Barrier(comm)
    msg = (; _sum=sum(f, A), _length=length(A))
    (; _sum, _length) = MPI.Allreduce(msg, _mean, comm)
    _sum / _length
end

@inline function _mean(a, b)
    (; _sum=(a._sum + b._sum), _length=(a._length + b._length))
end

"""
    allreduce(A, op, topo::Topology)

Element-wise reduction of `A` across all processes using binary operator `op`.
"""
allreduce(A, op, topo::Topology) = _allreduce(processing(topo), op, A)

@inline _allreduce(::SingleProcessing, op, A) = A

@inline function _allreduce(::MultiProcessing, op, A)
    comm = Topologies.comm()
    MPI.Barrier(comm)
    MPI.Allreduce(A, op, comm)
end

"""
    allgather(x, topo::Topology)

Gather value from all processes into a vector.
"""
function allgather(x, topo::Topology)
    if Topologies.is_multi_processing(topo)
        comm = Topologies.comm()
        MPI.Allgather(x, comm)
    else
        [x]
    end
end

"""
    barrier(topo::Topology)

Synchronize all processes at a barrier point.
"""
function barrier(topo::Topology)
    if Topologies.is_multi_processing(topo)
        comm = Topologies.comm()
        MPI.Barrier(comm)
    end
end
