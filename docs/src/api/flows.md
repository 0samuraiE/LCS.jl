# Flows API

```@meta
CurrentModule = LCS.Flows
```

Incompressible flow field computation with forcing and coupling.

## Configuration

```@docs
FlowConfig
FlowParams
```

## Buffer and State

```@docs
FlowBuffer
FlowState
```

## Operators

```@docs
init!
rhs!
update!
couple!
forcing!
timestep
timestep!
```

## Initial Conditions

```@docs
Init
RandomFlow
IdealFlow
```

## Forcing Schemes

```@docs
Forcing
RCForcing
EnergyPreserveRCF
NoForcing
```

## Coupling Schemes

```@docs
Couple
HSMAC
```

## Statistics

```@docs
Stat
FlowStat
stat
Q
```
