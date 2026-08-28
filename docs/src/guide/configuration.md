# Configuration & Environment

## Environment Variables

LCS.jl recognizes the following environment variables:

| Variable | Purpose | Default | How to Set | Notes |
|----------|---------|---------|-----------|-------|
| `LCS_LOG_LEVEL` | Control simulation logging verbosity | `INFO` | `LCS_LOG_LEVEL=QUIET julia script.jl` | User-facing |
| `LCS_RESUME` | Resume from checkpoint when output directory exists | `false` | `LCS_RESUME=1 julia script.jl` | User-facing |
| `LCS_TEST_SKIP_MPI` | Skip MPI tests during test execution | `false` | `LCS_TEST_SKIP_MPI=1 julia test/runtests.jl` | Development only |
| `JULIA_REFERENCETESTS_UPDATE` | Update reference test baselines | `false` | `JULIA_REFERENCETESTS_UPDATE=true julia test/...` | Development only |

`LCS_RESUME`, `LCS_TEST_SKIP_MPI`, and `JULIA_REFERENCETESTS_UPDATE` are boolean flags: any truthy value (`1`, `true`, `yes`) enables the option.

**User-facing variables:**
- `LCS_LOG_LEVEL`: Accepts `INFO` (default), `QUIET` (suppress all logging), or `PROFILE` (enable fine-grained timing output). Useful for batch runs (`QUIET`) or performance analysis (`PROFILE`).
- `LCS_RESUME`: Automatically continue from the latest checkpoint instead of overwriting output

**Development variables:**
- `LCS_TEST_SKIP_MPI`: Skip MPI tests when running the test suite locally
- `JULIA_REFERENCETESTS_UPDATE`: Regenerate reference test baselines after intentional output changes

## Configuration Files

LCS.jl simulations are configured using YAML files with the `.lcs-yaml` extension. Configuration is loaded via [`LCSIO.readconfig`](@ref):

```julia
using LCS
config = LCSIO.readconfig("simulation.lcs-yaml")
```

### Simulation Modes

Two simulation modes are available, specified via the `mode.kind` field:
- `flow`: Flow-only simulation
- `flow-particle`: Coupled flow and Lagrangian particle simulation

## Configuration Schema

The full configuration schema is defined in YAML meta-schema files located in `lib/schema/schemas/`:

- `common.meta.yaml` - Backend, topology, grid, timestep, simulation
- `flow.meta.yaml` - Flow parameters, initialization, coupling, forcing
- `particle.meta.yaml` - Particle population, drag models, initialization

### Meta-Schema Notation

LCS.jl uses a custom meta-schema notation that compiles to JSON Schema Draft 7. This provides a concise YAML syntax for defining configuration structure and validation rules.

#### Primitive Types

Declare types using type keywords:

```yaml
name: string
count: integer
value: number
enabled: boolean
```

#### Fixed-Length Arrays

Append array size in parentheses for fixed-length arrays:

```yaml
dims: integer(3)        # 3-element integer array (grid dimensions)
proc_dims: integer(3)   # 3-element integer array (MPI topology)
bounds: number(2)       # 2-element number array
```

#### Objects

Dictionaries declare objects. All fields are required by default, and no additional properties are allowed:

```yaml
params:
  U0_m_s: number
  L0_m: number
  NU_m2_s: number
```

#### String Constants

Non-type strings declare string constants:

```yaml
backend:
  oneOf:
    - kind: cpu      # "cpu" is a string constant
    - kind: cuda     # "cuda" is a string constant
```

#### Combinators

**`oneOf` - Exclusive choice:** Exactly one variant must match. Use `kind` field for discrimination:

```yaml
init:
  oneOf:
    - kind: random-flow
      params:
        seed: integer
    - kind: ideal-flow
    - kind: restart
      params:
        file: string
```

**`allOf` - Intersection:** All schemas must be satisfied.

**`anyOf` - Union:** At least one schema must be satisfied.

**`not` - Negation:** Schema must not match.

### Generating Schema

The meta-schemas compile to JSON Schema via [`Schema.generate`](@ref):

```julia
using Schema
Schema.generate()  # Produces lcs-schema.json
```

## See Also

- [`LCSIO.readconfig`](@ref) - Configuration loading
- [`Schema.generate`](@ref) - Schema generation and meta-schema notation
