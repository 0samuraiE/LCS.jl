# Parallel

```@meta
CurrentModule = Parallel
```

The `Parallel` module provides backend-agnostic utilities for parallel execution on CPUs and GPUs.

## Parallel iteration with `foraxes`

`Parallel.foraxes` executes an N-dimensional loop over the specified index ranges using a `KernelAbstractions.jl` backend. The same loop body can be used for CPU and GPU execution by changing only the backend.

```julia
Parallel.foraxes(
    backend,
    (2:N-1, 2:N-1, 2:N-1),
) do i, j, k
    dUdt[i, j, k] =
        -(U[i+1, j, k] - U[i-1, j, k]) / (2 * dx)
end
```

Here, `backend` specifies the execution backend, such as a CPU or GPU backend, and the tuple `(2:N-1, 2:N-1, 2:N-1)` specifies the iteration range in each dimension.

On a CPU backend, the loop is executed using Julia threads when multiple threads are available. On a GPU backend, `foraxes` generates and launches a GPU kernel through `KernelAbstractions.jl`.

An optional `coloring` keyword can be used for colored iteration schemes. For example, red–black iteration can be specified as follows:

```julia
Parallel.foraxes(
    backend,
    axes;
    coloring=Parallel.RedBlack(),
) do i, j, k
    # computation
end
```

## API reference

```@docs
Parallel.foraxes
Parallel.NoColoring
Parallel.RedBlack
Parallel.RedBlackFast
Parallel.RedBlackBlock
Parallel.hostaccess
```
