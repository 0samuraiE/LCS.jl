using BenchmarkTools
using CUDA
using Parallel

const SUITE = BenchmarkGroup()

backends = Tuple{String,KA.Backend}[("cpu", KA.CPU())]
CUDA.has_cuda() && push!(backends, ("cuda", CUDA.CUDABackend()))

for (name, backend) in backends
    let group = addgroup!(SUITE, name)
        for (name, coloring) in [
            ("no", Parallel.NoColoring()),
            ("red-black", Parallel.RedBlack()),
            ("red-black-fast", Parallel.RedBlackFast()),
            ("red-black-block", Parallel.RedBlackBlock(2)),
        ]
            let group = addgroup!(group, name)
                A = KA.zeros(backend, Float64, 256, 256, 256)
                x = sin(1.0)
                y = cos(1.0)

                group["foraxes"] = @benchmarkable begin
                    Parallel.foraxes($backend, axes($A); coloring=$coloring) do i, j, k
                        @inbounds ($A)[i, j, k] = $x / $y
                    end
                    KA.synchronize($backend)
                end
            end
        end
    end
end
