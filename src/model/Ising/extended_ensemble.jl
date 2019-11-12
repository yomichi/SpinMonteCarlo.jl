function init_extended_ensemble(model::Ising, param::Parameter)
    dosobsname = param["Observable for Extended Ensemble"]
    if dosobsname == "Number of Parallel Bonds"
        return DoS(0, numbonds(model), 1)
    elseif dosobsname == "Energy"
        delta = param["Bin Size"]
        return DoS(-2.0, 2.0, delta)
    elseif dosobsname == "Magnetization"
        delta = param["Bin Size"]
        return DoS(-1.0, 1.0, delta)
    else
        error()
    end
end
