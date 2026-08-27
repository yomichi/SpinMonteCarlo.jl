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

## Observables
```@docs
MCObservable
ScalarObservable
VectorObservable
MCObservableSet
makeMCObservable!
SimpleObservable
SimpleVectorObservable
SimpleObservableSet
SimpleVectorObservableSet
binning
Jackknife
JackknifeVector
JackknifeSet
JackknifeVectorSet
jackknife
extrapolate_tau
extrapolate_stderror
```

## Snapshot

`runMC` writes spin configurations when `param["Snapshot Interval"]` is a
positive number of measurement MCS (thermalization steps are never sampled).
The file is `"$(param["Snapshot Filename Prefix"])_$(param["ID"]).txt"`, one
configuration per line, values separated by spaces. Nothing identifies the
model or lattice in the file, so keeping track of which run produced which file
is up to you.

The flattening follows Julia's column-major order, so the length of a row is
not always the number of sites:

| Model | Row length | Order | Values |
|---|---|---|---|
| `Ising` | `numsites` | site | `±1` |
| `Potts`, `Clock` | `numsites` | site | `1` to `Q` |
| `XY` | `numsites` | site | `σ ∈ [0,1)`, angle `θ = 2πσ` |
| `AshkinTeller` | `2 * numsites` | `σ₁ τ₁ σ₂ τ₂ …` | `±1` |

Two runs sharing an `"ID"` and prefix write to the same file, exactly as they
would share a checkpoint file.

```@docs
snapshot
save_snapshot
load_snapshots
```

## Removed API
```@docs
gen_snapshot!
gensave_snapshot!
load_snapshot
```

## Utility
```@docs
Parameter
convert_parameter
convert_parameter(::Ising, ::Parameter)
convert_parameter(::Potts, ::Parameter)
convert_parameter(::Clock, ::Parameter)
convert_parameter(::XY, ::Parameter)
convert_parameter(::AshkinTeller, ::Parameter)
convert_parameter(::QuantumXXZ, ::Parameter)
```
