module ParallelExt
using Parallel
using PolySerde

@variant Parallel.Coloring
@composite Parallel.NoColoring
@composite Parallel.RedBlack
@composite Parallel.RedBlackFast
@composite Parallel.RedBlackBlock
@kind Parallel.Coloring "no" Parallel.NoColoring
@kind Parallel.Coloring "red-black" Parallel.RedBlack
@kind Parallel.Coloring "red-black-fast" Parallel.RedBlackFast
@kind Parallel.Coloring "red-black-block" Parallel.RedBlackBlock
end
