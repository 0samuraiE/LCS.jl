# LCS.jl

**Lagrangian Cloud Simulator in Julia** — A single-source, multi-platform,
multiphase turbulence simulation model with MPI and GPU support.

[![CI](https://github.com/0samuraiE/LCS.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/0samuraiE/LCS.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://0samuraiE.github.io/LCS.jl/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19515891.svg)](https://doi.org/10.5281/zenodo.19515891)

## Overview

LCS.jl is a direct numerical simulation (DNS) model for turbulent particle-laden flows, implemented in [Julia](https://julialang.org/) and [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
A single source code runs on CPUs, NVIDIA GPUs, and AMD GPUs. The design is intended to support other backends (Apple Metal) through KernelAbstractions.jl, though these have not been tested.

For details on the model formulation, validation, and performance, see the paper (link to be added upon arXiv release).

## Installation Requirements

| Component | Tested version |
|-----------|----------------|
| Julia | 1.12.x |
| OpenMPI | 5.0.x |
| HDF5 | 1.14.x |
| CUDA (required for NVIDIA GPU) | 12.8.x |
| cuDNN (optional) | 9.8.x |

Other versions may work but are untested.

Tested platforms:

| Backend | Hardware |
|---------|----------|
| CPU | AMD EPYC 9654, Intel Xeon Gold 6326 |
| NVIDIA GPU | H100 SXM5, RTX A6000 |
| AMD GPU | AMD Radeon RX 9070 XT |

## Installation

### With Devcontainer (recommended)

The recommended environment is provided as a devcontainer. Open this repository in VS Code with the Dev Containers extension. Dependencies are installed automatically.

### Without a dev container

**Requirements:** Julia 1.12, MPI, HDF5, and optionally CUDA or ROCm

Julia, MPI, and HDF5 must be installed manually. Configure `LocalPreferences.toml` for system MPI and HDF5 — ensure that MPI is built with CUDA- or ROCm-aware and parallel HDF5 support as needed. Build scripts used in the dev container are available in [`.devcontainer/scripts`](.devcontainer/scripts) as a reference. Then run:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Configuration

Simulation parameters are specified in YAML configuration files. Examples are provided in the `cfgs/` directory.
A JSON schema is provided for editor support (autocompletion and validation).
Key parameters include grid size, particle count, time step, and backend selection. See the [Documentation](https://0samuraiE.github.io/LCS.jl/) for the full configuration reference.

To regenerate the schema, run:

```bash
julia --project -e 'using Schema; Schema.generate()'
```

## Running

Choose a backend through one of the self-contained entry-point environments:

| Backend | Environment | Entry point |
|---------|-------------|-------------|
| CPU | `LCSCPU` | `LCSCPU/simulate.jl` |
| NVIDIA GPU | `LCSCUDA` | `LCSCUDA/simulate.jl` |
| AMD GPU | `LCSAMDGPU` | `LCSAMDGPU/simulate.jl` |

Instantiate the selected environment once before the first run, for example:

```bash
julia --project=LCSCUDA -e 'using Pkg; Pkg.instantiate()'
```

Every entry point accepts a YAML configuration path and an optional Julia expression
that patches the parsed configuration before the simulation starts.

The default forcing method (`rc-energy`) is the Reduced Communication Forcing (Onishi et al., 2011), which uses a box mean filter implemented via the DNN library (cuDNN for NVIDIA, MIOpen for AMD). If the DNN library is unavailable or unsupported on your platform, use the constant-power linear forcing (`linear-power`; Lundgren, 2003; Rosales & Meneveau, 2005) instead, which does not require a DNN library and runs on any backend:

```yaml
forcing:
  kind: "linear-power"
  params:
    power: 0.5
```

### CPU (single process, with profiling)

```bash
LCS_LOG_LEVEL=PROFILE julia --project=LCSCPU LCSCPU/simulate.jl cfgs/N128.lcs-yaml
```

### NVIDIA GPU (single process)

```bash
LCS_LOG_LEVEL=INFO julia --project=LCSCUDA LCSCUDA/simulate.jl cfgs/N128.lcs-yaml
```

### AMD GPU (single process)

```bash
LCS_LOG_LEVEL=INFO julia --project=LCSAMDGPU LCSAMDGPU/simulate.jl cfgs/N128.lcs-yaml
```

### Multi-GPU (MPI)

```bash
mpirun --allow-run-as-root -n 4 \
    -x LD_LIBRARY_PATH \
    -x UCX_WARN_UNUSED_ENV_VARS="n" \
    -x UCX_ERROR_SIGNALS="SIGILL,SIGBUS,SIGFPE" \
    -x LCS_LOG_LEVEL=INFO \
    julia --project=LCSCUDA LCSCUDA/simulate.jl cfgs/N1500.lcs-yaml
```

## Testing

```bash
JULIA_CPU_THREADS=2 LCS_TEST_SKIP_MPI=0 julia --project -e 'using Pkg; Pkg.test(; test_args=["--quickfail"])'
```

## Reproducing Figures

Figures in the paper can be reproduced with the following commands. Input data are stored in the `data/` directory.

```bash
julia --project=plots -e "using Pkg; Pkg.instantiate()"
julia --project=plots ./plots/runall.jl
```

## License

MIT License © Taketo Tominaga, Ryo Onishi
