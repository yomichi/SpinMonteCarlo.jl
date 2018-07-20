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
postproc(::Clock,::Parameter,::MCObservableSet)
postproc(::QuantumXXZ,::Parameter,::MCObservableSet)
```

## Lattice

```@autodocs
Modules = [SpinMonteCarlo]
Pages = ["src/lattice.jl"]
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
