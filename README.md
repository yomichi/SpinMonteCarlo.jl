# SpinMonteCarlo.jl
Markov chain Monte Carlo solver for finite temperature problem of lattie spin system implemented by [Julia](https://julialang.org) language.

[Online manual](https://yomichi.github.io/SpinMonteCarlo.jl/latest)

# Install

``` julia
Pkg> add SpinMonteCarlo
```

# Simple example

[The following program](example/simple.jl) calculates temperature v.s. specific heat of the ferromagnetic Ising model on a $16\times 16$ square lattice by Swendsen-Wang algorithm.

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
    @printf("%f %.15f %.15f\n",
            T, mean(result["Specific Heat"]), stderror(result["Specific Heat"]))
end
```

# Implemented 

## Model
- Classical spin model
    - `Ising` model
    - `Q` state `Potts` model
        - order parameter defined as $M = n_1(Q-1)/Q  - (1-n_1)/Q$, where $n_1$ is the number density of $q=1$ spins.
    - `XY` model
    - `Q` state `Clock` model
    - `AshkinTeller` model
- Quantum spin model
    - spin-`S` `QuantumXXZ` model
        - $\mathcal{H} = \sum_{ij} [ Jz_{ij} S_i^z S_j^z + \frac{Jxy_{ij}}{2} (S_i^+ S_j^- + S_i^-S_j^+) ] - \sum_i \Gamma_i S_i^x$

## Lattice
- `chain lattice`
    - `L`
- `bond-alternating chain lattice`
    - `L`
- `square lattice`
    - `L * W`
- `J1J2 square lattice`
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
        - $\braket{m} := \braket{ M_\text{total}/N_\text{site} }$
    - `|Magnetization|`
        - $\braket{|m|} := \braket{|M_\text{total}/N_\text{site}|}$ 
    - `Magnetization^2`
        - $\braket{m^2} := \braket{(M_\text{total}/N_\text{site})^2}$
    - `Magnetization^4`
        - $\braket{m^4} := \braket{(M_\text{total}/N_\text{site})^4 }$
    - `Binder Ratio`
        - $U_{4,2} := \braket{m^4}/\braket{m^2}^2$
    - `Susceptibility`
        - $\chi := \partial_h \braket{m} = (N/T)(\braket{m^2} - \braket{m}^2)$
    - `Connected Susceptibility`
        - $\chi := (N/T)(\braket{m^2} - \braket{|m|}^2)$
    - `Energy`
        - $E := \braket{\mathcal{H}} = \braket{E_\text{total}}/N_\text{site}$
    - `Energy^2`
        - $E^2 := \braket{\mathcal{H}^2}$
    - `Specific Heat`
        - $C := \partial_\beta \braket{\mathcal{H}} = (N/T^2)(\braket{\mathcal{H}^2} - \braket{\mathcal{H}}^2)$
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
        - $\braket{m} := \braket{\sum_i S_i^z } / N_\text{site}$
    - `Magnetization^2`
        - $\braket{m^2}:= \braket{(\sum_i S_i^z)^2 } / N_\text{site}^2$
    - `Magnetization^4`
        - $\braket{m^4}:= \braket{(\sum_i S_i^z)^4 } / N_\text{site}^4$
    - `Binder Ratio`
        - $U_{4,2} := \braket{m^4}/\braket{m^2}^2$
    - `Susceptibility`
        - $\chi := \partial_h \braket{m} = (N/T)(\braket{m^2} - \braket{m}^2)$
    - `Energy`
        - $E := \braket{\mathcal{H}}$
    - `Energy^2`
        - $E^2 := \braket{\mathcal{H}^2}$
    - `Specific Heat`
        - $C := \partial_\beta \braket{\mathcal{H}} = (N/T^2)(\braket{\mathcal{H}^2} - \braket{\mathcal{H}}^2)$

# Future work
- `Model`
    - Classical model
        - Heisenberg model
        - antiferro interaction
        - magnetic field
    - Quantum model
        - SU(N) model
- `Update Method`
    - worm algorithm
- Others
    - random number parallelization
        - NOTE: parameter parallelization can be realized simply by using `@parallel for` or `pmap`.
    - write algorithmic note
        - especially, Foutuin-Kasteleyn representaion and improved estimators

# Author
[Yuichi Motoyama](https://github.com/yomichi), the University of Tokyo, 2016-

This package is distributed under the MIT license.
