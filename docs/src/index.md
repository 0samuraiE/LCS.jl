# LCS.jl

Lagrangian Cloud Simulator in Julia with MPI and GPU support.

```julia
using MPI
MPI.Init()

using LCS
LCS.simulate("config.lcs-yaml")

MPI.Finalize()
```

## Guides

- **[Configuration & Environment](guide/configuration.md)**: Environment variables and configuration files

## Core Modules

- **[Flows](api/flows.md)**: Flow field computation
- **[Particles](api/particles.md)**: Lagrangian particle transport
- **[LCSIO](api/lcsio.md)**: Configuration and I/O
