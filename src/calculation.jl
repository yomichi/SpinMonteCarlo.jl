function runMC(params::Dict)
    verbose = get(params, "Verbose", false)
    if verbose
        println("Start: ", params)
    end
    model = params["Model"](params)
    T = params["T"]
    MCS = get(params, "MCS", 8192)
    Therm = get(params, "Thermalization", MCS>>3)
    blocal = get(params, "LocalUpdate", false)
    ret = runMC(model, T, MCS, Therm, blocal)
    if verbose
        println("Finish: ", params)
    end
    return ret
end

function runMC(model::Union{Ising, Potts}, T::Real, MCS::Integer, Therm::Integer, blocal::Bool)
    if blocal
        for mcs in 1:Therm
            local_update!(model,T)
        end
    else
        for mcs in 1:Therm
            SW_update!(model,T)
        end
    end

    nsites = numsites(model.lat)
    invV = 1.0/nsites
    obs = BinningObservableSet()
    makeMCObservable!(obs, "Magnetization")
    makeMCObservable!(obs, "|Magnetization|")
    makeMCObservable!(obs, "Magnetization^2")
    makeMCObservable!(obs, "Magnetization^4")
    makeMCObservable!(obs, "Energy")
    makeMCObservable!(obs, "Energy^2")

    if blocal
        for mcs in 1:MCS
            local_update!(model,T)
            M, E = measure(model, T)
            M2 = M*M
            M4 = M2*M2
            E *= invV
            obs["Magnetization"] << M
            obs["|Magnetization|"] << abs(M)
            obs["Magnetization^2"] << M2
            obs["Magnetization^4"] << M4
            obs["Energy"] << E
            obs["Energy^2"] << E*E
        end
    else
        for mcs in 1:MCS
            sw_info = SW_update!(model,T)
            M, M2, M4 = magnetizations(sw_info, model)
            E, E2 = energy(sw_info, model, T)
            obs["Magnetization"] << 0.0
            obs["|Magnetization|"] << abs(M)
            obs["Magnetization^2"] << M2
            obs["Magnetization^4"] << M4
            obs["Energy"] << E
            obs["Energy^2"] << E2
        end
    end

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites/T)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites/T)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites/T/T)*(jk["Energy^2"] - jk["Energy"]^2)

    return jk
end

function runMC(model::Union{Clock, XY}, T::Real, MCS::Integer, Therm::Integer, blocal::Bool)
    if blocal
        for i in 1:Therm
            local_update!(model, T)
        end
    else
        for i in 1:Therm
            SW_update!(model, T)
        end
    end

    obs = BinningObservableSet()
    makeMCObservable!(obs, "|Magnetization|")
    makeMCObservable!(obs, "|Magnetization|^2")
    makeMCObservable!(obs, "|Magnetization|^4")
    makeMCObservable!(obs, "|Magnetization x|")
    makeMCObservable!(obs, "Magnetization x^2")
    makeMCObservable!(obs, "Magnetization x^4")
    makeMCObservable!(obs, "|Magnetization y|")
    makeMCObservable!(obs, "Magnetization y^2")
    makeMCObservable!(obs, "Magnetization y^4")
    makeMCObservable!(obs, "Helicity Modulus x")
    makeMCObservable!(obs, "Helicity Modulus y")
    makeMCObservable!(obs, "Energy")
    makeMCObservable!(obs, "Energy^2")

    nsites = numsites(model.lat)
    invV = 1.0/nsites
    beta = 1.0/T

    if blocal
        for i in 1:MCS
            local_update!(model, T)
            measure_impl!(obs, model, T, invV)
        end
    else
        for i in 1:MCS
            SW_update!(model, T)
            measure_impl!(obs, model, T, invV)
        end
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

    return jk
end

function measure_impl!(obs, model::Union{Clock, XY}, T, invV)
    M, E, U = measure(model, T)
    E *= invV
    x2 = M[1]*M[1]
    y2 = M[2]*M[2]
    m2 = x2+y2
    x4 = x2*x2
    y4 = y2*y2
    m4 = m2*m2
    obs["|Magnetization x|"] << abs(M[1])
    obs["Magnetization x^2"] << x2
    obs["Magnetization x^4"] << x4
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

function print_result(io::IO, params, obs)
    i  = 1
    for name in ["Q", "L", "T"]
        println(io, "# $i : $name")
        i += 1
    end
    ks = collect(keys(obs[1]))
    sort!(ks)
    for name in ks
        println(io, "# $i, $(i+1) : $name")
        i += 2
    end

    for (p, o) in zip(params, obs)
        print(io, get(p, "Q", 0), " ")
        print(io, p["L"], " ")
        print(io, p["T"], " ")

        for name in ks
            print(io, mean(o[name]), " ")
            print(io, stderror(o[name]), " ")
        end
        println(io)
    end
end

