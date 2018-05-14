using JLD2

"""
    runMC(param::Dict)
    runMC(params::AbstractArray{T}; parallel::Bool=false) where T<:Dict

Runs Monte Carlo simulation and returns calculated observables.

If `"\$(param["Checkpoint Filename Prefix"])__\$(param["ID"]).jld2"` exists and
`param["Checkpoint Interval"] > 0.0`, this loads the checkpoint file and restarts the pending simulation.

When `parallel==true`, `runMC(params)` uses `pmap` instead of `map`.

# Required keys in `param`
- "Model"
    - `param` will be used as the argument of the constructor.
- "Update Method"
    - `param` will be used as an argument.
- "T": Temperature
# Optional keys in `param`
- "MCS": The number of Monte Carlo steps after thermalization
    - Default: 8192
- "Thermalization": The number of Monte Carlo steps for thermalization
    - Default: `MCS>>3`
- "Seed": The initial seed of the random number generator, `MersenneTwister`
    - Default: determined randomly (see `Random.srand`)
- "Checkpoint Filename Prefix"
    - Default: "cp"
- "Checkpoint Interval": Interval (in seconds) between saving checkpoint files
    - Default: 0.0, this means that NO checkpoint files will be loaded and saved.

"""
function runMC(params::AbstractArray{T}; parallel::Bool=false) where T<:Dict
    map_fn = ifelse(parallel, pmap, map)
    return map_fn(enumerate(params)) do id, p
        p["ID"] = id
        return runMC(p)
    end
end

function runMC(param::Parameter)
    model = param["Model"](param)
    if "Seed" in keys(param)
        srand(model, param["Seed"])
    end
    ret = runMC(model, param)
    return ret
end

function runMC(model, param::Parameter)
    verbose = get(param, "Verbose", false)
    if verbose
        println("Start: ", param)
    end
    cp_filename = @sprintf("%s_%d.jld2", get(param, "Checkpoint Filename Prefix", "cp"), get(param, "ID", 1))
    cp_interval = get(param, "Checkpoint Interval", 0.0)
    tm = time()

    MCS = get(param, "MCS", 8192)
    Therm = get(param, "Thermalization", MCS>>3)

    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")

    if cp_interval > 0.0 && ispath(cp_filename)
        @load(cp_filename, model, obs, mcs)
    end

    if "UpdateMethod" in keys(param)
        warn("\"UpdateMethod\" is deprecated. Use instead \"Update Method\".")
        param["Update Method"] = param["UpdateMethod"]
    end
    update! = param["Update Method"]
    estimator = get(param, "Estimator", default_estimator(model, update!))

    while mcs < MCS
        t = @elapsed begin
            st = update!(model,param)
            localobs = estimator(model, param, st)
        end
        if mcs >= Therm
            obs["Time per MCS"] << t
            accumulateObservables!(model, obs, localobs)
        end
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end

    if cp_interval > 0.0
        @save(cp_filename, model, obs, mcs)
    end

    jk = postproc(model, param, obs)

    if verbose
        println("Finish: ", param)
    end
    return jk
end

function accumulateObservables!(::Model, obs::MCObservableSet, localobs::Measurement)
    for key in keys(localobs)
        if haskey(obs, key)
            obs[key] << localobs[key]
        else
            makeMCObservable!(obs, key)
            obs[key] << localobs[key]
        end
    end
    return obs
end

function postproc(model::Union{Ising, Potts}, param, obs)
    nsites = numsites(model)
    T = param["T"]
    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0/jk["Time per MCS"]
    return jk
end

function postproc(model::Union{Clock, XY}, param, obs)
    nsites = numsites(model)
    T = param["T"]
    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio x"] = jk["Magnetization x^4"] / (jk["Magnetization x^2"]^2)
    jk["Binder Ratio y"] = jk["Magnetization y^4"] / (jk["Magnetization y^2"]^2)
    jk["Binder Ratio"] = jk["|Magnetization|^4"] / (jk["|Magnetization|^2"]^2)
    jk["Susceptibility x"] = (nsites*beta)*jk["Magnetization x^2"]
    jk["Susceptibility y"] = (nsites*beta)*jk["Magnetization y^2"]
    jk["Susceptibility"] = (nsites*beta)*jk["|Magnetization|^2"]
    jk["Connected Susceptibility x"] = (nsites*beta)*(jk["Magnetization x^2"] - jk["|Magnetization x|"]^2)
    jk["Connected Susceptibility y"] = (nsites*beta)*(jk["Magnetization y^2"] - jk["|Magnetization y|"]^2)
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["|Magnetization|^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0 / jk["Time per MCS"]
    return jk
end

function postproc(model::QuantumXXZ, param, obs)
    nsites = numsites(model)
    T = param["T"]
    beta = 1.0/T

    jk = jackknife(obs)

    for oname in ("Magnetization", "|Magnetization|",
                  "Magnetization^2", "Magnetization^4",
                  "Energy", "Energy^2",
                 )
        jk[oname] = jk["Sign * $oname"] / jk["Sign"]
    end

    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0/jk["Time per MCS"]
    return jk
end


