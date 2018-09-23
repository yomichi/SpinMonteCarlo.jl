# using JLD2

#=
@doc """
    runMC(param::Parameter)
    runMC(params::AbstractArray{Parameter}
          ;
          parallel::Bool=false,
          autoID::Bool=true)

Runs Monte Carlo simulation(s) and returns calculated observables.

If a checkpoint file named `"\$(param["Checkpoint Filename Prefix"])_\$(param["ID"]).jld2"` exists and
`param["Checkpoint Interval"] > 0.0`, `runMC` loads this file and restarts the pending simulation.

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
    - Default: determined randomly (see `Random.srand`)
- "Checkpoint Filename Prefix": See above document.
    - Default: `"cp"`
- "ID": See above document.
    - Default: `0`
- "Checkpoint Interval": Time interval between writing checkpoint file in seconds.
    - Default: `0.0`, this means that NO checkpoint file will be loaded and saved.
"""
=#
@doc """
    runMC(param::Parameter)
    runMC(params::AbstractArray{Parameter} ; parallel::Bool=false)

Runs Monte Carlo simulation(s) and returns calculated observables.

# Keyward aruguments
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
    - Default: determined randomly (see `Random.srand`)
"""
function runMC(params::AbstractArray{T}; parallel::Bool=false, autoID::Bool=true) where T<:Dict
    map_fn = ifelse(parallel, pmap, map)
    # return map_fn(enumerate(params)) do idp
    return map_fn(params) do p
        #=
        id,p = idp
        if autoID
            p["ID"] = id
        end
        =#
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
    #=
    cp_filename = @sprintf("%s_%d.jld2",
                           get(param, "Checkpoint Filename Prefix", "cp")::String,
                           get(param, "ID", 0)::Int)
    cp_interval = get(param, "Checkpoint Interval", 0.0) :: Float64
    =#
    tm = time()

    MCS = get(param, "MCS", 8192) :: Int
    Therm = get(param, "Thermalization", MCS>>3) :: Int

    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "MCS per Second")

    #=
    if cp_interval > 0.0 && ispath(cp_filename)
        jld = jldopen(cp_filename)
        model = jld["model"] :: MODEL
        obs = jld["obs"] :: BinningObservableSet
        mcs = convert(Int, jld["mcs"]) :: Int
        close(jld)
    end
    =#

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
            obs["MCS per Second"] << (1.0/t)
            accumulateObservables!(model, obs, localobs)
        end
        mcs += 1
        #=
        if cp_interval > 0.0 && time() - tm > cp_interval
            jld = jldopen(cp_filename, "w")
            jld["model"] = model
            jld["obs"] = obs
            jld["mcs"] = mcs
            close(jld)
            tm += cp_interval
         end
         =#
    end

    #=
    if cp_interval > 0.0
        jld = jldopen(cp_filename, "w")
        jld["model"] = model
        jld["obs"] = obs
        jld["mcs"] = mcs
        close(jld)
    end
    =#

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

@doc """
    postproc(model::Model, param::Parameter, obs::MCObservableSet)

Post process of observables. For example, Specific heat will be calculated from energy, energy^2, and temperature.
"""
function postproc end

@doc doc"""
    postproc(model::Union{Ising, Potts}, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $m$ is total magnetization per site and $\epsilon$ is total energy per site.

- `"Binder Ratio"`
    - $R := \frac{\left \langle m^4 \right \rangle}{\left \langle m^2 \right\rangle^2}$
- `"Susceptibility"`
    - $\chi := \frac{N}{T}\left(\left\langle m^2\right\rangle\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::Union{Ising, Potts}, param::Parameter, obs::MCObservableSet)
    nsites = numsites(model)
    T = param["T"] :: Float64
    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    return jk
end

@doc doc"""
    postproc(model::Union{Clock, XY}, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $m$ is total magnetization per site and $\epsilon$ is total energy per site.

- `"Binder Ratio x"`
    - $\frac{\left \langle m_x^4 \right \rangle}{\left \langle m_x^2 \right\rangle^2}$
- `"Binder Ratio y"`
    - $\frac{\left \langle m_y^4 \right \rangle}{\left \langle m_y^2 \right\rangle^2}$
- `"Binder Ratio"`
    - $\frac{\left \langle |m|^4 \right \rangle}{\left \langle |m|^2 \right\rangle^2}$
- `"Susceptibility x"`
    - $\frac{N}{T}\left(\left\langle m_x^2\right\rangle\right)$
- `"Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle m_y^2\right\rangle\right)$
- `"Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle |m|^2\right\rangle\right)$
- `"Connected Susceptibility x"`
    - $\frac{N}{T}\left(\left\langle m_x^2\right\rangle - \left\langle |m_x| \right\rangle^2\right)$
- `"Connected Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle m_y^2\right\rangle - \left\langle |m_y| \right\rangle^2\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle |m|^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::Union{Clock, XY}, param::Parameter, obs::MCObservableSet)
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
    return jk
end

@doc doc"""
    postproc(model::QuantumXXZ, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $s$ is sign of weight, $m$ is total magnetization per site,
and $\epsilon$ is total energy per site.

- `"Magnetization"`
    - $\left\langle m s\right\rangle\Big/\left\langle s \right\rangle$
- `"|Magnetization|"`
    - $\left\langle |m| s\right\rangle\Big/\left\langle s \right\rangle$
- `"Magnetization^2"`
    - $\left\langle m^2 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Magnetization^4"`
    - $\left\langle m^4 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Energy"`
    - $\left\langle \epsilon s\right\rangle\Big/\left\langle s \right\rangle$
- `"Energy^2"`
    - $\left\langle \epsilon^2 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Binder Ratio"`
    - $\frac{\left \langle m^4 \right \rangle}{\left \langle m^2 \right\rangle^2}$
- `"Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::QuantumXXZ, param::Parameter, obs::MCObservableSet)
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
    return jk
end
