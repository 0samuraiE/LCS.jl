module LCS
export LCSIO
export Flows
export Particles

using PolySerde
using MPI
using Offsets
using Parallel
using PPrint
using Printf
using Topologies
using Topologies.Utils

const N_DIMS = 3
const DOMAIN_ORIGIN = 0.0
const DOMAIN_LENGTH = 2π
const EMPTY_LOG = NamedTuple()

const FP = Float64
const IP = Int64

const ENV_LCS_RESUME = "LCS_RESUME"

const ENV_LCS_LOG_LEVEL = "LCS_LOG_LEVEL"
const LCS_LOG_QUIET = "QUIET"
const LCS_LOG_INFO = "INFO"
const LCS_LOG_PROFILE = "PROFILE"

get_log_level() = get(ENV, ENV_LCS_LOG_LEVEL, LCS_LOG_INFO)

include("aliases.jl")
using .Aliases

include("utils/Utils.jl")
using .Utils
using .Utils: @ntuple, @concretize

include("rungekutta.jl")
include("config.jl")
include("buffer.jl")
include("state.jl")
include("diff.jl")

include("topologies.jl")

include("log.jl")
include("init.jl")
include("evolve.jl")
include("timestep.jl")
include("stat.jl")
include("simulate.jl")

include("flows/Flows.jl")
include("particles/Particles.jl")
include("LCSIO/LCSIO.jl")

include("precompile.jl")
end
