const Parameter = Dict{String, Any}

function convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)
    T = param["T"]
    J = get(param, "J", 1.0)
    if isa(J, Real)
        Js = J.*ones(numbondtypes(model))
    else
        Js = Vector{Float64}(J)
    end
    return T, Js
end

function convert_parameter(model::QuantumXXZ, param::Parameter)
    T = param["T"]

    Jz = get(param, "Jz", 1.0)
    if isa(Jz, Real)
        Jzs = Jz.*ones(numbondtypes(model))
    else
        Jzs = Vector{Float64}(Jz)
    end

    Jxy = get(param, "Jxy", 1.0)
    if isa(Jxy, Real)
        Jxys = Jxy.*ones(numbondtypes(model))
    else
        Jxys = Vector{Float64}(Jxy)
    end

    G = get(param, "Gamma", 0.0)
    if isa(G, Real)
        Gs = G.*ones(numsitetypes(model))
    else
        Gs = Vector{Float64}(G)
    end

    return T, Jzs, Jxys, Gs
end

