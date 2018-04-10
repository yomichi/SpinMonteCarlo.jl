using JLD2

function runMC(params::AbstractArray; parallel::Bool=false)
    map_fn = ifelse(parallel, pmap, map)
    return map_fn(enumerate(params)) do id, p
        p["ID"] = id
        return runMC(p)
    end
end

function runMC(param::Dict)
    verbose = get(param, "Verbose", false)
    cp_filename = @sprintf("%s_%d.jld2", get(param, "Checkpoint Filename Prefix", "cp"), get(param, "ID", 1))
    if verbose
        println("Start: ", param)
    end
    model = param["Model"](param)
    if "Seed" in keys(param)
        srand(model, param["Seed"])
    end
    ret = runMC(model, param, cp_filename=cp_filename, cp_interval=get(param, "Checkpoint Interval", 0.0))
    if verbose
        println("Finish: ", param)
    end
    return ret
end

function runMC(model::Union{Ising, Potts}, param::Dict; cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    T = param["T"]
    Js = param["J"]
    MCS = get(param, "MCS", 8192)
    Therm = get(param, "Thermalization", MCS>>3)
    update! = get(param, "UpdateMethod", SW_update!)
    return runMC(model, T, Js, MCS, Therm, update!, cp_filename=cp_filename, cp_interval=cp_interval)
end
function runMC(model::Union{Ising, Potts}, T::Real, Js::Union{Real,AbstractArray}, MCS::Integer, Therm::Integer, update! = SW_update!
               ;
               cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    tm = time()
    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "Magnetization")
    makeMCObservable!(obs, "|Magnetization|")
    makeMCObservable!(obs, "Magnetization^2")
    makeMCObservable!(obs, "Magnetization^4")
    makeMCObservable!(obs, "Energy")
    makeMCObservable!(obs, "Energy^2")


    if cp_interval > 0.0 && ispath(cp_filename)
        @load(cp_filename, model, obs, mcs)
    end


    while mcs < Therm
        update!(model,T,Js,measure=false)
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
            localobs = update!(model,T,Js)
        end
        M = localobs["M"]
        M2 = localobs["M2"]
        M4 = localobs["M4"]
        E = localobs["E"]
        E2 = localobs["E2"]

        obs["Time per MCS"] << t
        obs["Magnetization"] << M
        obs["|Magnetization|"] << abs(M)
        obs["Magnetization^2"] << M2
        obs["Magnetization^4"] << M4
        obs["Energy"] << E
        obs["Energy^2"] << E2
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end

    if cp_interval > 0.0
        @save(cp_filename, model, obs, mcs)
    end

    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0/jk["Time per MCS"]

    return jk
end

function runMC(model::Union{Clock, XY}, param::Dict; cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    T = param["T"]
    Js = param["J"]
    MCS = get(param, "MCS", 8192)
    Therm = get(param, "Thermalization", MCS>>3)
    update! = get(param, "UpdateMethod", SW_update!)
    return runMC(model, T, Js, MCS, Therm, update!, cp_filename=cp_filename, cp_interval=cp_interval)
end
function runMC(model::Union{Clock, XY}, T::Real, Js::Union{Real,AbstractArray}, MCS::Integer, Therm::Integer, update! =SW_update!
               ;
               cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    tm = time()
    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
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

    if cp_interval > 0.0 && ispath(cp_filename)
        @load(cp_filename, model, obs, mcs)
    end

    while mcs < Therm
        update!(model,T,Js,measure=false)
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end


    nsites = numsites(model)
    invV = 1.0/nsites
    beta = 1.0/T

    while mcs < MCS
        t = @elapsed begin
            localobs = update!(model, T, Js)
        end

        M = localobs["M"]
        E = localobs["E"]
        U = localobs["U"]

        x2 = M[1]*M[1]
        y2 = M[2]*M[2]
        m2 = x2+y2
        x4 = x2*x2
        y4 = y2*y2
        m4 = m2*m2
        obs["Time per MCS"] << t
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
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end
    if cp_interval > 0.0
        @save(cp_filename, model, obs, mcs)
    end

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

function runMC(model::QuantumXXZ, param::Dict
               ;
               cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    T = param["T"]
    if "J" in keys(param)
        Jz = Jxy = param["J"]
    else
        Jz = param["Jz"]
        Jxy = param["Jxy"]
    end
    G = get(param, "Gamma", 0.0)
    MCS = get(param, "MCS", 8192)
    Therm = get(param, "Thermalization", MCS>>3)
    return runMC(model, T, Jz, Jxy, G, MCS, Therm, cp_filename=cp_filename, cp_interval=cp_interval)
end
function runMC(model::QuantumXXZ, T::Real,
               Jz::Union{Real, AbstractArray}, Jxy::Union{Real, AbstractArray},
               Gs::Union{Real, AbstractArray}, MCS::Integer, Therm::Integer
               ;
               cp_filename::AbstractString="cp.jld2", cp_interval::Real=0.0)
    tm = time()
    mcs = 0
    MCS += Therm
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "Sign * Magnetization")
    makeMCObservable!(obs, "Sign * |Magnetization|")
    makeMCObservable!(obs, "Sign * Magnetization^2")
    makeMCObservable!(obs, "Sign * Magnetization^4")
    makeMCObservable!(obs, "Sign * Energy")
    makeMCObservable!(obs, "Sign * Energy^2")
    makeMCObservable!(obs, "Sign")

    if cp_interval > 0.0 && ispath(cp_filename)
        @load(cp_filename, model, obs, mcs)
    end

    while mcs < Therm
        loop_update!(model,T, Jz, Jxy, Gs, measure=false)
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end

    nsites = numsites(model.lat)
    invV = 1.0/nsites
    while mcs < MCS
        t = @elapsed begin 
            localobs = loop_update!(model,T,Jz,Jxy,Gs)
        end
        M = localobs["M"]
        M2 = localobs["M2"]
        M4 = localobs["M4"]
        E = localobs["E"]
        E2 = localobs["E2"]
        sgn = localobs["Sign"]
        obs["Time per MCS"] << t
        obs["Sign * Magnetization"] << M*sgn
        obs["Sign * |Magnetization|"] << abs(M)*sgn
        obs["Sign * Magnetization^2"] << M2*sgn
        obs["Sign * Magnetization^4"] << M4*sgn
        obs["Sign * Energy"] << E*sgn
        obs["Sign * Energy^2"] << E2*sgn
        obs["Sign"] << sgn
        mcs += 1
        if cp_interval > 0.0 && time() - tm > cp_interval
            @save(cp_filename, model, obs, mcs)
            tm += cp_interval
        end
    end
    if cp_interval > 0.0
        @save(cp_filename, model, obs, mcs)
    end

    jk = jackknife(obs)

    for oname in ("Magnetization", "|Magnetization|",
                  "Magnetization^2", "Magnetization^4",
                  "Energy", "Energy^2",
                 )
        jk[oname] = jk["Sign * $oname"] / jk["Sign"]
    end

    beta = 1.0/T

    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    jk["MCS per Second"] = 1.0/jk["Time per MCS"]

    return jk
end

