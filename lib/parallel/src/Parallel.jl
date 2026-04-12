module Parallel
export KA, foraxes

using KernelAbstractions

import KernelAbstractions as KA

include("foraxes.jl")
include("hostaccess.jl")
end
