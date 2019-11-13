function init_extended_ensemble(model::Ising, param::Parameter)
    dosobsname = param["Observable for Extended Ensemble"]
    if dosobsname == "Number of Parallel Bonds"
        return DoS(0, numbonds(model), 1)
    elseif dosobsname == "Energy"
        E = 0.0
        Js = convert_parameter(model, param)[2]
        for bt in 1:numbondtypes(model)
            J = abs(Js[bt])
            nbonds = numbonds(model, bt)
            E += nbonds * J
        end
        E /= numsites(model)

        delta = param["Bin Size"]
        return DoS(-E, E, delta)
    elseif dosobsname == "Magnetization"
        delta = param["Bin Size"]
        return DoS(-1.0, 1.0, delta)
    else
        error()
    end
end
