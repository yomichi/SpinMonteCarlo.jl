export learnEE, measureEE

function updateEE!(model::Model, present, dos::DoS, dosobsname, p)
    action = nextaction(model)
    diff = localchange(model, action, p...)[dosobsname]
    if rand() < exp(logg(dos, present) - logg(dos, present+diff))
        present += diff
        accept!(model, action)
    end
    return present
end


function learnEE(param::Parameter)
    model = param["Model"](param)
    return learnEE(model, param)
end

function learnEE(model::Model, param::Parameter)
    dos = init_extended_ensemble(model, param)
    dosobsname = param["Observable for Extended Ensemble"]

    p = convert_parameter(model, param)

    α = get(param, "WL Start Refinement Factor", 1.0)
    α_last = get(param, "WL Minimum Refinement Factor", 1e-4)

    criteria = get(param, "WL Criteria for Flat", 0.8)
    nonzero_ratio = get(param, "WL Nonzero Ratio", 0.4)

    interval_calc_full_state = 100

    istage = 1

    while α >= α_last
        reset_hist!(dos)
        iter = 0
        present = simple_estimator(model, p..., nothing)[dosobsname]
        while true
            for i in 1:numsites(model)
                present = updateEE!(model, present, dos, dosobsname, p)
                visit!(dos, present, α)
            end
            if check_flat(dos, α, criteria = criteria, nonzero_ratio = nonzero_ratio, verbose=true)
                break
            end
            iter += 1
            if iter % interval_calc_full_state == 0
                present = simple_estimator(model, p..., nothing)[dosobsname]
                iter = 0
            end
        end
        renormalize!(dos)
        open("dos_$(istage).dat", "w") do io
            if valuetype(dos) <: Integer
                for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                    v = index2value(dos, i)
                    @printf(io, "%d %.15f %d\n", v, g, h)
                end
            else
                for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                    @printf(io, "%.15f %.15f %d\n", v, g, h)
                end
            end
        end
        α *= 0.5
        istage += 1
    end

    save_dos("dos.jld", dos)

    return dos
end

function measureEE(model::Model, dos::DoS, param::Parameter)
    MCS = get(param, "MCS", 65536)
    Thermalization = get(param, "Thermalization", MCS >> 3)

    dosobsname = param["Observable for Extended Ensemble"]
    p = convert_parameter(model, param)
    present = simple_estimator(model, p..., nothing)[dosobsname]

    for mcs in 1:Thermalization
        for i in 1:numsites(model)
            present = updateEE!(model, present, dos, dosobsname, p)
        end
        present = simple_estimator(model, p..., nothing)[dosobsname]
    end

    obs = BinningObservableSet()

    for mcs in 1:MCS
        for i in 1:numsites(model)
            present = updateEE!(model, present, dos, dosobsname, p)
        end
        lobs = simple_estimator(model, p..., nothing)
        lobs["Identity"] = 1.0
        present = lobs[dosobsname]
        accumulateObservablesEE!(model, logg(dos, present), obs, lobs)
    end
    
    jk = jackknife(obs)
    for key in keys(jk)
        if ! startswith(key, "Wg ")
            continue
        end
        obsname = string(key[4:end])
        jk[obsname] = jk[key] / jk["Wg Identity"]
    end

    jk = postproc(model, param, jk)

    return jk
end

@doc """
    accumulateObservablesEE!(model, logg, obs::MCObservableSet, localobs::Dict)

Accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["Energy"]`.
"""
function accumulateObservablesEE!(::Model, logg, obs::MCObservableSet, localobs::Measurement)
    if length(obs) < 3
        @inbounds for key in keys(localobs)
            obsname = "Wg " * key
            makeMCObservable!(obs, obsname)
            v = localobs[key]
            v *= exp(localobs["Log Boltzmann Weight"] + logg)
            obs[obsname] << v
        end
    else
        @inbounds for key in keys(localobs)
            v = localobs[key]
            v *= exp(localobs["Log Boltzmann Weight"] + logg)
            obs["Wg " * key] << v
        end
    end
    return obs
end
