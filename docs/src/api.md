# API

## Driver

```@meta
CurrentModule = SpinMonteCarlo
```

```@docs
runMC
initObservables
accumulateObservables!
postproc
```

## Model

```@docs
Ising
Potts
Clock
XY
QuantumXXZ
```

## Lattice

```@autodocs
Modules = [SpinMonteCarlo]
Pages = ["src/lattice.jl"]
```

## Lattice generator

```@docs
dimer_lattice
chain_lattice
square_lattice
triangular_lattice
cubic_lattice
fully_connected_lattice
```

