# SpinMonteCarlo.jl
[![Build Status](https://travis-ci.org/yomichi/SpinMonteCarlo.jl.svg)](https://travis-ci.org/yomichi/SpinMonteCarlo.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ooxw9wusg26bklq3?svg=true)](https://ci.appveyor.com/project/yomichi/spinmontecarlo-jl)

Markov chain Monte Carlo solver for finite temperature problem of lattie spin system implemented by [Julia](https://julialang.org) language.

# Install

``` julia
julia> Pkg.clone("https://github.com/yomichi/SpinMonteCarlo.jl")
```

# Simple example

[The following program](example/ising.jl) calculates temperature v.s. specific heat of the ferromagnet Ising model on a 16x16 square lattice by Swendsen-Wang algorithm.

``` julia
using SpinMonteCarlo

const model = Ising
const lat = square_lattice
const L = 16
const update = SW_update!

const Tc = 2.0/log1p(sqrt(2))
const Ts = Tc*linspace(0.85, 1.15, 31)
const MCS = 8192
const Therm = MCS >> 3

for T in Ts
    params = Dict{String,Any}( "Model"=>model, "Lattice"=>lat,
                                 "L"=>L, "T"=>T,
                                 "UpdateMethod"=>update,
                                 "MCS"=>MCS, "Thermalization"=>Therm,
                             )
    result = runMC(params)
    println(@sprintf("%f %.15f %.15f",
                      T, mean(result["Specific Heat"]), stderror(result["Specific Heat"])))
end
```

# Implemented 

## Model
- Classical spin model
    - `Ising` model
    - `Q` state `Potts` model
        - order parameter defined as `M = (Q-1)/Q * n_1 - (1-n_1)/Q`, where `n_1` is the number density of `q=1` spins and `N` is the number of all spins.
    - `XY` model
    - `Q` state `Clock` model
- Quantum spin model
    - S=1/2 `TransverseFieldIsing` model

## Lattice
- `chain_lattice`
- `square_lattice`
- `triangular_lattice`
- `cubic_lattice`

## Update algorithm
- Classical spin
    - `local_update!`
    - `SW_update!`
    - `Wolff_update!`
- Quantum spin
    - `loop_update!`

## Physical quantities
- `Ising`, `Potts`
    - `Magnetization`
    - `|Magnetization|`
    - `Magnetization^2`
    - `Magnetization^4`
    - `Binder Ratio`
    - `Susceptibility`
    - `Connected Susceptibility`
    - `Energy`
    - `Energy^2`
    - `Specific Heat`
- `XY`, `Clock`
    - `|Magnetization|`
    - `|Magnetization|^2`
    - `|Magnetization|^4`
    - `Binder Ratio`
    - `Susceptibility`
    - `Connected Susceptibility`
    - `Magnetization x`
    - `|Magnetization x|`
    - `Magnetization x^2`
    - `Magnetization x^4`
    - `Binder Ratio x`
    - `Susceptibility x`
    - `Connected Susceptibility x`
    - `Magnetization y`
    - `|Magnetization y|`
    - `Magnetization y^2`
    - `Magnetization y^4`
    - `Binder Ratio y`
    - `Susceptibility y`
    - `Connected Susceptibility y`
    - `Helicity Modulus x`
    - `Helicity Modulus y`
    - `Energy`
    - `Energy^2`
    - `Specific Heat`
- `TransverseFieldIsing`
    - `Magnetization`
    - `|Magnetization|`
    - `Magnetization^2`
    - `Magnetization^4`
    - `Binder Ratio`
    - `Susceptibility`
    - `Connected Susceptibility`

# Future work
- `Model`
    - Classical model
        - Heisenberg model
        - antiferro interaction
        - magnetic field
    - Quantum model
        - Heisenberg model
        - XXZ model
        - SU(N) model
- `Lattice`
    - ladder
    - tube
    - anisotropic parameters (e.g. J1-J2)
- `UpdateMethod`
    - worm algorithm
- Others
    - resume and restart
    - random number parallelization
        - NOTE: parameter parallelization can be realized simply by using `@parallel for` or `pmap`.
    - write algorithmic note
        - especially, Foutuin-Kasteleyn representaion and improved estimators

# Author
[Yuichi Motoyama](https://github.com/yomichi), the University of Tokyo

This package distributed under the MIT license.
