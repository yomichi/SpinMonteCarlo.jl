# Internal APIs

Documentation for `SpinMonteCarlo.jl`'s internals (not exported).

## Driver

```@meta
CurrentModule = SpinMonteCarlo
```

```@docs
accumulateObservables!
postproc
postproc(::Ising,::Parameter,::MCObservableSet)
postproc(::Potts,::Parameter,::MCObservableSet)
postproc(::Clock,::Parameter,::MCObservableSet)
postproc(::XY,::Parameter,::MCObservableSet)
postproc(::AshkinTeller,::Parameter,::MCObservableSet)
postproc(::QuantumXXZ,::Parameter,::MCObservableSet)
```

## Lattice

```@docs
generatelattice
numsitetypes
numbondtypes
```

## Model

```@docs
LoopElementType
LET_Cut
LET_FMLink
LET_AFLink
LET_Vertex
LET_Cross
LocalLoopOperator
```

## Utility
```@docs
default_estimator
@gen_convert_parameter
SWInfo
UnionFind
addnode!
unify!
clusterize!
clusterid
root!
root_path_halving!
root_path_splitting!
```
