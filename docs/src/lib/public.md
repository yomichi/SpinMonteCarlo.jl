# Public interfaces

Documentation for `SpinMonteCarlo.jl`'s public interface (exported).

## Driver

```@meta
CurrentModule = SpinMonteCarlo
```

```@docs
runMC
```

## Model

```@docs
Ising
Potts
Clock
XY
QuantumXXZ
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

## Update method

An index of model parameter (e.g., `Js`) is corresponding to `sitetype` or `bondtype`.

```@docs
local_update!
SW_update!
Wolff_update!
loop_update!
```

## Estimator
```@docs
simple_estimator
improved_estimator
```

## Utility
```@docs
convert_parameter
convert_parameter(::Ising, ::Parameter)
convert_parameter(::QuantumXXZ, ::Parameter)
```
