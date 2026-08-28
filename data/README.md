# Plot data

This directory contains the normalized input data used by the plotting scripts
in [`plots/`](../plots). From the repository root, regenerate every figure with:

```bash
julia --project=plots -e 'using Pkg; Pkg.instantiate()'
julia --project=plots plots/runall.jl
```

## Contents

| Files | Description |
|-------|-------------|
| `cpu-strong.csv`, `cpu-weak.csv` | CPU strong- and weak-scaling results |
| `cpu-strong-intra.csv` | Intra-node CPU strong-scaling results |
| `gpu-strong.csv`, `gpu-weak.csv` | GPU strong- and weak-scaling results |
| `comm-optimization.csv` | Halo-exchange optimization results |
| `gpu-cpu-compare.csv` | GPU/CPU performance comparison |
| `julia-vs-fortran.csv` | Julia/Fortran performance comparison |
| `hetero.csv` | Heterogeneous execution results |
| `flow-stats.csv` | Turbulence statistics |
| `energy-spectrum-*.h5` | Energy-spectrum statistics |
| `rdf-at-contact.csv` | Reference and present RDF-at-contact values |
| `rdf-model-parameters.csv` | RDF reference-model parameters |

The Julia/Fortran comparison uses an N=384 case for both implementations. The
Julia values cover 1–8 MPI ranks; the Fortran values use 100 warm-up steps and
100 measured steps.

GPU weak-scaling efficiency uses the 8-GPU result as the common baseline for
all device counts. The GPU scaling summary is authoritative where measurements
from the same configuration overlap.
