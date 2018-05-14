const Parameter = Dict{String, Any}

function convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)
    T = Float64(param["T"])
    nbt = numbondtypes(model)

    J = get(param, "J", 1.0) :: Union{Float64, AbstractArray{Float64}}
    Js = zeros(nbt)
    Js .= J
    return T, Js
end

function convert_parameter(model::QuantumXXZ, param::Parameter)
    T = Float64(param["T"])

    nbt = numbondtypes(model)
    nst = numsitetypes(model)

    Jz = get(param, "Jz", 1.0)
    Jzs = zeros(nbt)
    Jzs .= Jz

    Jxy = get(param, "Jxy", 1.0)
    Jxys = zeros(nbt)
    Jxys .= Jxy

    G = get(param, "Gamma", 0.0)
    Gs = zeros(nbt)
    Gs .= G

    return T, Jzs, Jxys, Gs
end

