# LCS Core API

```@meta
CurrentModule = LCS
```

Simulation entry point, configuration types, and time integration.

## Simulation

```@docs
simulate
```

## Configuration

```@docs
Config
Grid
TimeStep
SimulateConfig
TopoDraft
```

## State Management

```@docs
State
```

## Buffers

```@docs
AbstractBuffer
Buffer
```

## Time Integration

```@docs
timestep
evolve!
init!
```

## Statistics

```@docs
stat
printstate
printstat
printlog
```

## Modes

```@docs
Mode
FlowMode
FlowParticleMode
```
