module Utils
export mpiprintln, @mpishow

using MPI
using Printf
using Statistics
using Topologies

mpiprintln(msgs...) = mpiprintln(stdout, msgs...)
function mpiprintln(io::IO, msgs...)
    root = 0

    comm = Topologies.comm()
    rank = MPI.Comm_rank(comm)

    msgss = MPI.gather(msgs, comm; root)
    if rank == root
        for (i, msgs) in enumerate(msgss)
            println(io, "[rank $(i - 1)] ", msgs...)
        end
    end

    MPI.Barrier(comm)
end

macro mpishow(exs...)
    Expr(
        :block,
        [
            :(mpiprintln($(sprint(Base.show_unquoted, ex)), " = ", Base.repr(
                begin
                    local value = $(esc(ex))
                end,
            ))) for ex in exs
        ]...,
    )
end

struct Summary
    mean :: Float64
    min  :: Float64
    max  :: Float64
    std  :: Float64
end

function summary(elapsed::Real, comm::MPI.Comm)
    elapseds = MPI.Allgather(elapsed, comm)
    Summary(mean(elapseds), minimum(elapseds), maximum(elapseds), std(elapseds))
end

function summarize(suite)
    for (name, summaries) in suite
        vals = getfield.(summaries, :mean)
        m = mean(vals)
        s = std(vals)
        println("$name: $m $s")
    end
end

macro mpibench(n, block)
    quote
        comm = Topologies.comm()
        rank = MPI.Comm_rank(comm)
        nprocs = MPI.Comm_size(comm)

        summaries = Vector{Summary}(undef, $(esc(n)))
        for itr in 1:($(esc(n)))
            GC.gc()
            GC.gc()
            GC.gc()
            GC.enable(false)
            MPI.Barrier(comm)
            elapsed = MPI.Wtime()
            $(esc(block))
            elapsed = MPI.Wtime() - elapsed
            GC.enable(true)

            s = summary(elapsed, comm)
            if rank == 0
                @printf "Iteration %d:\n" itr
                @printf "ave : %9.6f seconds\n" s.mean
                @printf "min : %9.6f seconds\n" s.min
                @printf "max : %9.6f seconds\n" s.max
                @printf "std : %9.6f seconds\n" s.std
            end

            summaries[itr] = summary(elapsed, comm)
        end

        summaries
    end
end

macro root(ex)
    quote
        comm = Topologies.comm()
        rank = MPI.Comm_rank(comm)

        if rank == 0
            $(esc(ex))
        end
    end
end
end
