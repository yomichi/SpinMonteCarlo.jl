# Internal APIs

Documentation for `SpinMonteCarlo.jl`'s internals.

## Driver

```@meta
CurrentModule = SpinMonteCarlo
```

```@docs
accumulateObservables!
postproc
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
UnionFind
addnode!
unify!
clusterize!
clusterid
```
