# SpinMonteCarlo.jl
<!--[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://yomichi.github.io/SpinMonteCarlo.jl/stable)-->

[![Docs (latest)](https://img.shields.io/badge/docs-latest-blue.svg)](https://yomichi.github.io/SpinMonteCarlo.jl/latest)
[![Build status (macos/Linux)](https://travis-ci.org/yomichi/SpinMonteCarlo.jl.svg)](https://travis-ci.org/yomichi/SpinMonteCarlo.jl)
[![Build status (win)](https://ci.appveyor.com/api/projects/status/ooxw9wusg26bklq3?svg=true)](https://ci.appveyor.com/project/yomichi/spinmontecarlo-jl)

Markov chain Monte Carlo solver for finite temperature problem of lattie spin system implemented by [Julia](https://julialang.org) language.

# Install

``` julia
julia> Pkg.add("SpinMonteCarlo")
```

# Simple example

[The following program](example/ising.jl) calculates temperature v.s. specific heat of the ferromagnetic Ising model on a 16x16 square lattice by Swendsen-Wang algorithm.

``` julia
using SpinMonteCarlo
using Printf

const model = Ising
const lat = "square lattice"
const L = 16
const update = SW_update!

const Tc = 2.0/log1p(sqrt(2))
const Ts = Tc*range(0.85, stop=1.15, length=31)
const MCS = 8192
const Therm = MCS >> 3

for T in Ts
    params = Dict{String,Any}("Model"=>model, "Lattice"=>lat,
                              "L"=>L, "T"=>T, "J"=>1.0,
                              "Update Method"=>update,
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
        - order parameter defined as `M = (Q-1)/Q * n_1 - (1-n_1)/Q`, where `n_1` is the number density of `q=1` spins.
    - `XY` model
    - `Q` state `Clock` model
- Quantum spin model
    - spin-`S` `QuantumXXZ` model
        - `\mathcal{H} = \sum_{ij} [ Jz_{ij} S_i^z S_j^z + \frac{Jxy_{ij}}{2} (S_i^+ S_j^- + S_i^-S_j^+) ] - \sum_i Gamma_i S_i^x`

## Lattice
- `chain lattice`
    - `L`
- `bond-alternating chain lattice`
    - `L`
- `square lattice`
    - `L * W`
- `triangular lattice`
    - `L * W`
- `cubic lattice`
    - `L * W * H`
- `fully connected graph`
    - `N`

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
        - `:= < M_total/N_site >`
    - `|Magnetization|`
        - `:= < |M_total/N_site| >`
    - `Magnetization^2`
        - `:= < (M_total/N_site)^2 >`
    - `Magnetization^4`
        - `:= < (M_total/N_site)^4 >`
    - `Binder Ratio`
        - `:= "Magnetization^4")/"Magnetization^2"^2`
    - `Susceptibility`
        - `:= \partial_h "Magnetization" = (N/T) * ("Magnetization^2" - "Magnetization"^2)`
    - `Connected Susceptibility`
        - `:= (N_site/T) * ("Magnetization^2"- "|Magnetization|"^2)`
    - `Energy`
        - `:= < E_total/N_site >`
    - `Energy^2`
        - `:= < (E_total/N_site)^2 >`
    - `Specific Heat`
        - `:= \partial_beta "Energy" = (N/T^2)*("Energy^2" - "Energy"^2)`
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
- `QuantumXXZ`
    - `Magnetization`
        - `:= < \sum_i S_i^z > / Nsite`
    - `Magnetization^2`
        - `:= < (\sum_i S_i^z)^2 > / Nsite^2`
    - `Magnetization^4`
        - `:= < (\sum_i S_i^z)^4 > / Nsite^4`
    - `Binder Ratio`
        - `:= "Magnetization^4")/"Magnetization^2"^2`
    - `Susceptibility`
        - `:= \partial_h "Magnetization" = (N/T) * ("Magnetization^2" - "Magnetization"^2)`
    - `Energy`
        - `:= < \mathcal{H} > / Nsite`
    - `Energy^2`
        - `:= < \mathcal{H}^2 > / Nsite^2`
    - `Specific Heat`
        - `:= \partial_beta "Energy" = (N/T^2)*("Energy^2" - "Energy"^2)`

# Future work
- `Model`
    - Classical model
        - Heisenberg model
        - antiferro interaction
        - magnetic field
    - Quantum model
        - SU(N) model
- `Lattice`
    - open boundary condition
- `UpdateMethod`
    - worm algorithm
- Others
    - random number parallelization
        - NOTE: parameter parallelization can be realized simply by using `@parallel for` or `pmap`.
    - write algorithmic note
        - especially, Foutuin-Kasteleyn representaion and improved estimators

# Author
[Yuichi Motoyama](https://github.com/yomichi), the University of Tokyo, 2016-

This package distributed under the MIT license.
