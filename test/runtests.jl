using LCS
using Test
using ParallelTestRunner

const ENV_LCS_TEST_SKIP_MPI = "LCS_TEST_SKIP_MPI"

testsuite = Dict(
    "base" => quote
        include("test_base.jl")
    end,
    "utils" => quote
        include("test_utils.jl")
    end,
    "flows" => quote
        include("test_flows.jl")
    end,
    "particles" => quote
        include("test_particles.jl")
    end,
    "LCSIO" => quote
        include("test_LCSIO.jl")
    end,
    "integration" => quote
        include("test_integration.jl")
    end,
)

if !Base.get_bool_env(ENV_LCS_TEST_SKIP_MPI, false)
    push!(testsuite, "mpi" => quote
        include("test_mpi.jl")
    end)
end

ENV[LCS.ENV_LCS_LOG_LEVEL] = "QUIET"
runtests(LCS, ARGS; testsuite)
delete!(ENV, LCS.ENV_LCS_LOG_LEVEL)
