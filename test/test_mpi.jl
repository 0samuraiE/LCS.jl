using Test

using MPI

@testset "MPI" begin
    proc_dims = (3, 2, 2)

    expr = """
    using MPI
    using Test

    MPI.Init()

    include("$(joinpath(@__DIR__, "test_flows_mpi.jl"))")
    include("$(joinpath(@__DIR__, "test_particles_mpi.jl"))")
    include("$(joinpath(@__DIR__, "test_LCSIO_mpi.jl"))")
    include("$(joinpath(@__DIR__, "test_integration_mpi.jl"))")

    @test_nowarn MPI.Finalize()
    """

    cmd = `$(MPI.mpiexec()) --allow-run-as-root -n $(prod(proc_dims)) --map-by :oversubscribe $(Base.julia_cmd()) -e $expr -- $(join(proc_dims, " "))`
    p = run(ignorestatus(cmd))
    @test success(p)
end
