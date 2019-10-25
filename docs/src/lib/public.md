# Public interfaces

Documentation for `SpinMonteCarlo.jl`'s public interface (exported).

## Driver

```@meta
CurrentModule = SpinMonteCarlo
```

```@docs
runMC
```

## Lattice

```@docs
dim
size
sites
bonds
numsites
numbonds
neighbors
neighborsites
neighborbonds
source
target
sitetype
bondtype
sitecoordinate
bonddirection
cellcoordinate
```

## Model

```@docs
Ising
Potts
Clock
XY
AshkinTeller
QuantumXXZ
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
convert_parameter(::Potts, ::Parameter)
convert_parameter(::Clock, ::Parameter)
convert_parameter(::XY, ::Parameter)
convert_parameter(::AshkinTeller, ::Parameter)
convert_parameter(::QuantumXXZ, ::Parameter)
```
