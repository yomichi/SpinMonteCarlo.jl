# Run Monte Carlo

## Overview

`runMC` is a driver routine for a normal MCMC scheme having the following steps:

- Initialization
    - Initialize a model
- Thermalization
    - Update the spin configuration w/o measurement
- Measurement
    - Update the spin configuration w/  measurement

`runMC` takes a `param :: Parameter` and returns a MCMC result, `jks ::JackknifeSet` .
`param` is a set of parameters.
`jks` is a set of MCMC averaged quantities.

### Input
`param::Parameter` is a mapping from keynames to values (`Parameter` is an alias of `Dict{String, Any}`).
`runMC` requires the following parameters:

- "Model"
    - typename of a model (ex. `Ising`)
- "Lattice"
    - name of a lattice (ex. `"square lattice"`)
- "Update Method"
    - function to update a model (ex. `local_update!`)

The following are optional:

- "MCS"
    - The number of Monte Carlo steps for measurement
    - default: 8192
- "Thermalization"
    - The number of Monte Carlo steps for thermalization
    - default: "MCS" >> 3
- "Seed"
    - The initial seed of the random number generator, `Random.MersenneTwister`
    - default: random (see the doc of `Random.seed!`)
- "Checkpoint Filename Prefix"
    - Suffix of filename of checkpoint file (see the "Restart" section)
    - default: "cp"
- "ID"
    - Job ID used for restarting
    - default: 0
- "Checkpoint Interval"
    - Time interval between saving a calculation state into the checkpoint file in units of second.
    - default: 0.0, this means NO checkpoint files will be loaded and saved.
- "Verbose"
    - print message with the parameter before and after calculation
    - default: false

### Output
`jks::JackknifeSet` is a mapping from names to MCMC observables (`JackknifeSet` is an alias of `Dict{String, Jackknife}`).
`jk::Jackknife` is a MCMC observable whose statistics are calculated by [jack-knife method](https://en.wikipedia.org/wiki/Jackknife_resampling).
Users can retrieve mean, variance, standard deviation, standard error, and ``p \%`` confidence interval of `jk` by `mean(jk)`, `var(jk)`, `stddev(jk)`, `stderror(jk)`, and `confidence_interval(jk, p/100)`, respectively.

The set of calculated quantities depends on the model (an estimator and a postproc).
Please see the online document of specific functions (ex., `simple_estimator(::Ising)`).

### Multiple simulations
`runMC` can take a list of parameters as an input instead of a single parameter and performs MCMC calculation for each parameter and returns a list of results.
When a keyword arg `parallel = true` is passed, `runMC` performs MCMC calculations in parallel by using `pmap` function instead of `map`.

### Restart
If `param["Checkpoint Interval"] > 0.0`, `runMC` saves the state of calculation (the state of the model and the random number generator) into a checkpoint file named `"$(param["Checkpoint Filename Prefix"])_$(param["ID"]).dat"` every `param["Checkpoint Interval"]` seconds.

If a checkpoint file exists and `param["Checkpoint Interval"] > 0.0`, `runMC` loads this file and restarts the pending simulation.

`runMC(params::AbstractArray)` has a keyword argument, `autoID::Bool` (default=`true`).
If `true`, `["ID"]` will be automatically set (overwritten) as `params[i]["ID"] = i`.

NOTE: Restart will fail if the version or the system image of julia changes (see the doc of `Serialization.serialize`).

