using Serialization
using Distributed
import Distributed.pmap

export runMC

@doc """
    runMC(param::Parameter)
    runMC(params::AbstractArray{Parameter}
          ;
          parallel::Bool=false,
          autoID::Bool=true)

Runs Monte Carlo simulation(s) and returns calculated observables.

# Restart

If a checkpoint file named `"\$(param["Checkpoint Filename Prefix"])_\$(param["ID"]).dat"` exists and
`param["Checkpoint Interval"] > 0.0`, `runMC` loads this file and restarts the pending simulation.
NOTE: Restart will fail if the version or the system image of julia change (see the doc of `Serialization.serialize` ).

# Keyward aruguments
- `autoID`: If true, `"ID"`s will be set (overwritten) as `params[i]["ID"] = i`.
- `parallel`: If true, runs simulations in parallel (uses `pmap` instead of `map`).

# Required keys in `param`
- "Model"
- "Update Method"

# Optional keys in `param`
- "MCS": The number of Monte Carlo steps after thermalization
    - Default: `8192`
- "Thermalization": The number of Monte Carlo steps for thermalization
    - Default: `MCS>>3`
- "Seed": The initial seed of the random number generator, `MersenneTwister`
    - Default: determined randomly (see the doc of `Random.seed!`)
- "Checkpoint Filename Prefix": See the "Restart" section.
    - Default: `"cp"`
- "ID": See the "Restart" section.
    - Default: `0`
- "Checkpoint Interval": Time interval between writing checkpoint file in seconds.
    - Default: `0.0`, this means that NO checkpoint file will be loaded and saved.
"""
function runMC(params::AbstractArray{T}; parallel::Bool=false, autoID::Bool=true) where T<:Dict
    map_fn = ifelse(parallel, pmap, map)
    return map_fn(enumerate(params)) do (id,p)
        if autoID
            p["ID"] = id
        end
        return runMC(p)
    end
end

function runMC(param::Parameter)
    model = param["Model"](param)
    if "Seed" in keys(param)
        seed!(model, param["Seed"])
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
    cp_filename = @sprintf("%s_%d.dat",
                           get(param, "Checkpoint Filename Prefix", "cp")::String,
                           get(param, "ID", 0)::Int)
    cp_interval = get(param, "Checkpoint Interval", 0.0) :: Float64
    tm = time()

    MCS = get(param, "MCS", 8192) :: Int
    Therm = get(param, "Thermalization", MCS>>3) :: Int

    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "MCS per Second")

    if cp_interval > 0.0 && ispath(cp_filename)
        open(cp_filename) do io
            model = deserialize(io)
            obs = deserialize(io)
            mcs = deserialize(io)
        end
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
            obs["MCS per Second"] << (1.0/t)
            accumulateObservables!(model, obs, localobs)
        end
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            open(cp_filename, "w") do io
                serialize(io, model)
                serialize(io, obs)
                serialize(io, mcs)
            end
            tm += cp_interval
         end
    end

    if cp_interval > 0.0
        open(cp_filename, "w") do io
            serialize(io, model)
            serialize(io, obs)
            serialize(io, mcs)
        end
    end

    jk = postproc(model, param, obs)

    if verbose
        println("Finish: ", param)
    end
    return jk
end

@doc """
    accumulateObservables!(model, obs::MCObservableSet, localobs::Dict)

Accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["Energy"]`.
"""
function accumulateObservables!(::Model, obs::MCObservableSet, localobs::Measurement)
    if length(obs) < 3
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
