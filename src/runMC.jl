using JLD2

"""
    runMC(param::Dict)
    runMC(params::AbstractArray{T}; parallel::Bool=false) where T<:Dict

Runs Monte Carlo simulation and returns calculated observables.

If a checkpoint file named `"\$(param["Checkpoint Filename Prefix"])_\$(param["ID"]).jld2"` exists and
`param["Checkpoint Interval"] > 0.0`, `runMC` loads this file and restarts the pending simulation.

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
    MODEL = typeof(model)
    verbose = get(param, "Verbose", false) :: Bool
    if verbose
        println("Start: ", param)
    end
    cp_filename = @sprintf("%s_%d.jld2", get(param, "Checkpoint Filename Prefix", "cp")::String, get(param, "ID", 1)::Int)
    cp_interval = get(param, "Checkpoint Interval", 0.0) :: Float64
    tm = time()

    MCS = get(param, "MCS", 8192) :: Int
    Therm = get(param, "Thermalization", MCS>>3) :: Int

    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")

    if cp_interval > 0.0 && ispath(cp_filename)
        jld = jldopen(cp_filename)
        model = jld["model"] :: MODEL
        obs = jld["obs"] :: BinningObservableSet
        mcs = convert(Int, jld["mcs"]) :: Int
        close(jld)
    end

    if "UpdateMethod" in keys(param)
        warn("\"UpdateMethod\" is deprecated. Use instead \"Update Method\".")
        param["Update Method"] = param["UpdateMethod"] :: Function
    end
    update! = param["Update Method"] :: Function
    estimator = get(param, "Estimator", default_estimator(model, update!)) :: Function
    p = convert_parameter(model, param)

    while mcs < MCS
        if mcs < Therm
            update!(model, p...)
        else
            t = @elapsed begin
                st = update!(model,p...)
                localobs = estimator(model, p..., st)
            end
            obs["Time per MCS"] << t
            accumulateObservables!(model, obs, localobs)
        end
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            jld = jldopen(cp_filename, "w")
            jld["model"] = model
            jld["obs"] = obs
            jld["mcs"] = mcs
            close(jld)
            tm += cp_interval
         end
    end

    if cp_interval > 0.0
        jld = jldopen(cp_filename, "w")
        jld["model"] = model
        jld["obs"] = obs
        jld["mcs"] = mcs
        close(jld)
    end

    jk = postproc(model, param, obs)

    if verbose
        println("Finish: ", param)
    end
    return jk
end

"""
    accumulateObservables!(model, obs::MCObservableSet, localobs::Dict)

accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["Energy"]`.
"""
function accumulateObservables!(::Model, obs::MCObservableSet, localobs::Measurement)
    if length(obs) == 1
        @inbounds for key in keys(localobs)
            makeMCObservable!(obs, key)
            obs[key] << localobs[key]
        end
    else
        @inbounds for key in keys(localobs)
            obs[key] << localobs[key]
        end
    end
    return obs
end

"""
    postproc(model::Model, param::Dict, obs::MCObservableSet)

post process of observables. For example, Specific heat will be calculated from energy, energy^2, and temperature.
"""
function postproc end

function postproc(model::Union{Ising, Potts}, param, obs)
    nsites = numsites(model)
    T = param["T"] :: Float64
    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0/jk["Time per MCS"]
    return jk
end

function postproc(model::Union{Clock, XY}, param::Dict, obs::MCObservableSet)
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

function postproc(model::QuantumXXZ, param::Dict, obs::MCObservableSet)
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


