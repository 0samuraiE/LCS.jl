# Topologies

```@meta
CurrentModule = Topologies
```

MPI 3D Cartesian decomposition, halo exchange, and global reductions.

## Topology

```@docs
Topologies.Topology
Topologies.isroot
Topologies.device!
Topologies.proc_dims
Topologies.proc_size
Topologies.linear_rank
Topologies.linear_ranks
Topologies.cart_rank
Topologies.cart_to_linear_rank
Topologies.linear_to_cart_rank
Topologies.each_cart_rank
```

## Processing Mode

```@docs
Topologies.Processing
Topologies.SingleProcessing
Topologies.MultiProcessing
Topologies.processing
Topologies.is_multi_processing
```

## Halo Exchange

```@docs
Topologies.HaloBuffer
Topologies.SPHaloBuffer
Topologies.MPHaloBuffer
Topologies.HaloStyle
Topologies.FaceHalo
Topologies.FullHalo
Topologies.synchalo!
```

## Global Reduction

```@docs
Topologies.allsum
Topologies.allminimum
Topologies.allmaximum
Topologies.allmean
```
