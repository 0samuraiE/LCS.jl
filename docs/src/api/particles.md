# Particles API

```@meta
CurrentModule = LCS.Particles
```

Lagrangian particle advection with drag models and initialization.

## Configuration

```@docs
ParticleConfig
ParticleParams
Population
Cell
```

## Buffer and State

```@docs
ParticleBuffer
ParticleState
```

## Operators

```@docs
init!
motion!
update!
boundary!
timestep
timestep!
```

## Initialization

```@docs
Init
InitId
GenerateId
InitPosition
RandomPosition
InitVelocity
RestVelocity
InitSize
ConstSize
```

## Drag Models

```@docs
DragModel
LinearDrag
NonlinearDrag
```

## Statistics

```@docs
ParticleStat
stat
density
```
