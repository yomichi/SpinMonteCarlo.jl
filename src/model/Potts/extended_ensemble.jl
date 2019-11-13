function init_extended_ensemble(model::Potts, param::Parameter)
    dosobsname = param["Observable for Extended Ensemble"]
    if dosobsname == "Number of Parallel Bonds"
        return DoS(0, numbonds(model), 1)
    elseif dosobsname == "Energy"
        FME = 0.0
        AFE = 0.0
        Js = convert_parameter(model, param)[2]
        for bt in 1:numbondtypes(model)
            J = Js[bt]
            nbonds = numbonds(model, bt)
            if J > 0.0
                FME += nbonds * J
            else
                AFE -= nbonds * J
            end
        end
        FME /= numsites(model)
        AFE /= numsites(model)

        delta = param["Bin Size"]
        return DoS(-FME, AFE, delta)
    elseif dosobsname == "Magnetization"
        delta = param["Bin Size"]
        return DoS(-1.0/model.Q, 1.0-1.0/model.Q, delta)
    else
        error()
    end
end
