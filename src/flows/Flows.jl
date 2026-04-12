module Flows
using ..LCS
using ..LCS: Utils
using ..LCS.Aliases

using FFTW
using PolySerde
using MPI
using NNlib
using Offsets
using Parallel
using PPrint
using Random
using Statistics
using Topologies

import FourierTools as FT

include("config.jl")
include("buffer.jl")
include("state.jl")

include("inits/init.jl")
include("inits/idealflow.jl")
include("inits/randomflow.jl")

include("timestep.jl")
include("rhs.jl")
include("update.jl")

include("couples/couple.jl")
include("couples/hsmac.jl")

include("forcings/forcing.jl")
include("forcings/lf.jl")
include("forcings/rcf.jl")

include("stat.jl")
end
