# Develop Monte Carlo

## Lattice
See `src/lattice/standard.jl`.

## Model
`Model` should contain following fields: `lat :: Lattice` and `rng :: Random.MersenneTwister`.
`Model` also should have a constructor taking `param :: Parameter` as a argument.

You should define `convert_parameter` for your `Model`.
This is a helper function which takes `Model` and `Parameter`
and returns arguments of update method and estimator.

`@gen_convert_parameter` helps you to define `convert_parameter`.
For example, if your `model::A` needs a scalar `T = param["T"]` and a vector `Js = param["J"]` with `numbondtypes(model)` elements,

``` julia
@gen_convert_parameter(A, ("T", 1, 1.0), ("J", numbondtypes, 1.0))
```

defines your *documented and type-stable*  `convert_parameter(model::A, param::Parameter)` which returns `T` and `Js`.

### Note:

1. That the second element of each tuple is not a `Function` means that a return value is a scalar (the case of "T").
2. The third element of each tuple is the default value.

## Update method
"Update method" is a function which (in-place) updates `model::Model` under some parameters such as temperature `T`.
For example, `local_update!(model::Ising, T, Js)` updates a spin configuration of `model` by local spin flip and Metropolice-Hasting algorithm under temperature `T` and coupling constants `Js`.
"Update method" can return some object which will be used in "Estimator" as extra information.

### Example
Swendsen-Wang algorithm `SW_update!(::Ising, ::Parameter)` returns cluster information `sw::SWInfo` used in `improved_estimator`.

## Estimator
"Estimator" is a function which returns observables of a configuration `model` as a `Dict{String, Any}`.
Arguments of a "Estimator" are `model::Model`, parameters (return of `convert_parameter`), and extra information (return of "Update method") in order.

Default "Estimator" is determined by `model` and "Update Method" as return of `default_estimator(model, param["Update Method"])`.

### Example
`improved_estimator(model::Ising, T::Real, Js::AbstractArray, sw::SWInfo)` takes a return of `SW_update!(model::Ising, T::Real, Js::AbstractArray)` as the last argument.

## Postprocess
`postproc(model::Model, param::Parameter, obs::MCObservableSet)` is a post process of a Monte Carlo run.
Most important objective is to calculate functions of expectation values stored in `obs`.

For example, "Specific Heat" $C$ can be calculated from "Energy^2" $\langle E^2\rangle \big/ N^2$,
"Energy" $\left\langle E\right\rangle \big/ N$, the number of site $N$, and temperature $T$
as $C = \left(N\big/T^2\right)\left[ \left\langle E^2 \right\rangle - \left \langle E \right\rangle^2 \right]$.
This is realized as

``` julia
jk = jackknife(obs)
T = param["T"]
N = numsites(model)
jk["Specific Heat"] = (N/(T*T)) * (jk["Energy^2"] - jk["Energy"]^2)
```

`jackknife` converts `MCObservableSet` (e.g. `BinningObservableSet`) to `JackknifeObservableSet`,
which enables to calculate functions of mean values and these statistical errors.

`postproc` returns a `MCObservableSet` (usually `jk::JackknifeObservableSet` above),
which is also the return value of `runMC`.
