module LCSIO
using ..LCS
using ..LCS: Utils, @log
using ..LCS.Aliases

using Accessors
using HDF5
using PolySerde
using MPI
using Offsets
using OrderedCollections
using Parallel
using Printf
using Topologies
using Topologies.Utils
using YAML

const FILE_CONFIG = "config.lcs-yaml"
const FILE_STAT = "stat.h5"
const FILE_RESTART_REGEX = r"restart\.h5\.(\d+)$"
const FILE_RESTART_FMT = Printf.Format("restart.h5.%d")

const TAG_STATE = "state"
const TAG_STAT = "stat"
const TAG_LOG = "log"

const TAG_FLOW = "flow"
const TAG_FLOW_U = "U"
const TAG_FLOW_V = "V"
const TAG_FLOW_W = "W"
const TAG_FLOW_P = "P"

const TAG_PARTICLE = "particles"
const TAG_PARTICLE_COUNTS = "counts"
const TAG_PARTICLE_ID = "id"
const TAG_PARTICLE_DIAM = "diam"
const TAG_PARTICLE_X = "x"
const TAG_PARTICLE_Y = "y"
const TAG_PARTICLE_Z = "z"
const TAG_PARTICLE_U = "u"
const TAG_PARTICLE_V = "v"
const TAG_PARTICLE_W = "w"

include("utils.jl")
include("config.jl")
include("hdf5.jl")
include("save.jl")
include("load.jl")
include("upsample.jl")
end
