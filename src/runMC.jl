using JLD2

"""
    runMC(param::Dict)
    runMC(params::AbstractArray{T}; parallel::Bool=false) where T<:Dict

Runs Monte Carlo simulation and returns calculated observables.

If a file named `"\$(param["Checkpoint Filename Prefix"])_\$(param["ID"]).jld2"` exists and
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

function runMC(param::Dict)
    model = param["Model"](param)
    if "Seed" in keys(param)
        srand(model, param["Seed"])
    end
    ret = runMC(model, param)
    return ret
end

function runMC(model, param)
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
    obs = initObservables(model, binning=get(param,"Binning",true))

    if cp_interval > 0.0 && ispath(cp_filename)
        @load(cp_filename, model, obs, mcs)
    end

    if "UpdateMethod" in keys(param)
        warn("\"UpdateMethod\" is deprecated. Use instead \"Update Method\".")
        param["Update Method"] = param["UpdateMethod"]
    end
    update! = param["Update Method"]

    while mcs < Therm
        update!(model,param,measure=false)
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end

    nsites = numsites(model)
    invV = 1.0/nsites
    while mcs < MCS
        t = @elapsed begin
            localobs = update!(model,param)
        end
        obs["Time per MCS"] << t
        accumulateObservables!(model, obs, localobs)
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

"""
    initObservables(model::Model; binning::Bool=true)

    returns `BinningObservableSet` or `SimpleObservableSet` with some observables registered.
"""
function initObservables end

function initObservables(::Union{Ising, Potts}; binning::Bool=true)
    obs = ifelse(binning, BinningObservableSet(), SimpleObservableSet())
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "Magnetization")
    makeMCObservable!(obs, "|Magnetization|")
    makeMCObservable!(obs, "Magnetization^2")
    makeMCObservable!(obs, "Magnetization^4")
    makeMCObservable!(obs, "Energy")
    makeMCObservable!(obs, "Energy^2")
    return obs
end

function initObservables(::Union{Clock, XY}, binning::Bool=true)
    obs = ifelse(binning, BinningObservableSet(), SimpleObservableSet())
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "|Magnetization|")
    makeMCObservable!(obs, "|Magnetization|^2")
    makeMCObservable!(obs, "|Magnetization|^4")
    makeMCObservable!(obs, "Magnetization x")
    makeMCObservable!(obs, "|Magnetization x|")
    makeMCObservable!(obs, "Magnetization x^2")
    makeMCObservable!(obs, "Magnetization x^4")
    makeMCObservable!(obs, "Magnetization y")
    makeMCObservable!(obs, "|Magnetization y|")
    makeMCObservable!(obs, "Magnetization y^2")
    makeMCObservable!(obs, "Magnetization y^4")
    makeMCObservable!(obs, "Helicity Modulus x")
    makeMCObservable!(obs, "Helicity Modulus y")
    makeMCObservable!(obs, "Energy")
    makeMCObservable!(obs, "Energy^2")
    return obs
end

function initObservables(::QuantumXXZ; binning::Bool=true)
    obs = ifelse(binning, BinningObservableSet(), SimpleObservableSet())
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "Sign * Magnetization")
    makeMCObservable!(obs, "Sign * |Magnetization|")
    makeMCObservable!(obs, "Sign * Magnetization^2")
    makeMCObservable!(obs, "Sign * Magnetization^4")
    makeMCObservable!(obs, "Sign * Energy")
    makeMCObservable!(obs, "Sign * Energy^2")
    makeMCObservable!(obs, "Sign")
    return obs
end

"""
    accumulateObservables!(model, obs::MCObservableSet, localobs::Dict)

accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["E"]`.
"""
function accumulateObservables! end

function accumulateObservables!(model::Union{Ising,Potts}, obs::MCObservableSet, localobs::Dict)
    M = localobs["M"]
    M2 = localobs["M2"]
    M4 = localobs["M4"]
    E = localobs["E"]
    E2 = localobs["E2"]

    obs["Magnetization"] << M
    obs["|Magnetization|"] << abs(M)
    obs["Magnetization^2"] << M2
    obs["Magnetization^4"] << M4
    obs["Energy"] << E
    obs["Energy^2"] << E2
end

function accumulateObservables!(model::Union{Clock,XY}, obs::MCObservableSet, localobs::Dict)
    M = localobs["M"]
    E = localobs["E"]
    U = localobs["U"]

    x2 = M[1]*M[1]
    y2 = M[2]*M[2]
    m2 = x2+y2
    x4 = x2*x2
    y4 = y2*y2
    m4 = m2*m2
    obs["Magnetization x"] << M[1]
    obs["|Magnetization x|"] << abs(M[1])
    obs["Magnetization x^2"] << x2
    obs["Magnetization x^4"] << x4
    obs["Magnetization y"] << M[2]
    obs["|Magnetization y|"] << abs(M[2])
    obs["Magnetization y^2"] << y2
    obs["Magnetization y^4"] << y4
    obs["|Magnetization|"] << sqrt(m2)
    obs["|Magnetization|^2"] << m2
    obs["|Magnetization|^4"] << m4
    obs["Helicity Modulus x"] << U[1]
    obs["Helicity Modulus y"] << U[2]
    obs["Energy"] << E
    obs["Energy^2"] << E*E
end

function accumulateObservables!(model::QuantumXXZ, obs::MCObservableSet, localobs::Dict)
    M = localobs["M"]
    M2 = localobs["M2"]
    M4 = localobs["M4"]
    E = localobs["E"]
    E2 = localobs["E2"]
    sgn = localobs["Sign"]
    obs["Sign * Magnetization"] << M*sgn
    obs["Sign * |Magnetization|"] << abs(M)*sgn
    obs["Sign * Magnetization^2"] << M2*sgn
    obs["Sign * Magnetization^4"] << M4*sgn
    obs["Sign * Energy"] << E*sgn
    obs["Sign * Energy^2"] << E2*sgn
    obs["Sign"] << sgn
end

"""
    postproc(model::Model, param::Dict, obs::MCObservableSet)

post process of observables. For example, Specific heat will be calculated from energy, energy^2, and temperature.
"""
function postproc end

function postproc(model::Union{Ising, Potts}, param::Dict, obs::MCObservableSet)
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


