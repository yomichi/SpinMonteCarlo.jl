@doc doc"""
    default_estimator(model, updatemethod!)

Determines estimator to be used when `param["Estimator"]` is not set.
"""
default_estimator(model::Model, update) = simple_estimator
default_estimator(model::QuantumXXZ, update) = improved_estimator
default_estimator(model::Union{Ising, Potts}, update) = ifelse(update==SW_update!, improved_estimator, simple_estimator)

function simple_estimator(model, param::Parameter, _=nothing)
    p = convert_parameter(model,param)
    return simple_estimator(model, p..., nothing)
end

@doc """
    simple_estimator(model::Ising, T::Real, Js::AbstractArray)
    simple_estimator(model::Potts, T::Real, Js::AbstractArray)

Returns the following observables as `Dict{String, Any}`

# Observables
- `"Energy"`
    - Energy per spin (site)
- `"Energy^2"`
- `"Magnetization"`
    - Total magnetization per spin (order paremeter)
- `"|Magnetization|"`
- `"Magnetization^2"`
- `"Magnetization^4"`
"""
function simple_estimator(model::Ising, T::Real, Js::AbstractArray, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)

    M = mean(model.spins)
    E = 0.0
    @inbounds for b in bonds(model)
        s1, s2 = source(b), target(b)
        E += ifelse(model.spins[s1] == model.spins[s2], -1.0, 1.0) * Js[bondtype(b)]
    end
    E /= nsites

    res = Measurement()
    res["Magnetization"] = M
    res["|Magnetization|"] = abs(M)
    res["Magnetization^2"] = M^2
    res["Magnetization^4"] = M^4
    res["Energy"] = E
    res["Energy^2"] = E^2
    return res
end

function simple_estimator(model::Potts, T::Real, Js::AbstractArray, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)
    invQ = 1.0/model.Q

    M = 0.0
    @inbounds for s in 1:nsites
        M += ifelse(model.spins[s]==1, 1.0-invQ, -invQ)
    end
    M /= nsites
    E = 0.0
    @inbounds for b in bonds(model)
        s1, s2 = source(b), target(b)
        E -= ifelse(model.spins[s1] == model.spins[s2], 1.0, 0.0) * Js[bondtype(b)]
    end
    E /= nsites

    res = Measurement()
    res["Magnetization"] = M
    res["|Magnetization|"] = abs(M)
    res["Magnetization^2"] = M^2
    res["Magnetization^4"] = M^4
    res["Energy"] = E
    res["Energy^2"] = E^2
    return res
end


@doc """
    simple_estimator(model::Clock, T::Real, Js::AbstractArray)
    simple_estimator(model::XY, T::Real, Js::AbstractArray)

Returns the following observables as `Dict{String, Any}`

# Observables
- `"Energy"`
    - Energy per spin (site)
- `"Energy^2"`
- `"|Magnetization|"`
    - Absolute value of total magnetization per spin (order paremeter)
- `"|Magnetization|^2"`
- `"|Magnetization|^4"`
- `"Magnetization x"`
    - x component of total magnetization per spin (order paremeter)
- `"|Magnetization x|"`
- `"Magnetization x^2"`
- `"Magnetization x^4"`
- `"Magnetization y"`
    - y component of total magnetization per spin (order paremeter)
- `"|Magnetization y|"`
- `"Magnetization y^2"`
- `"Magnetization y^4"`
- `"Helicity Modulus x"`
- `"Helicity Modulus y"`
"""
function simple_estimator(model::Clock, T::Real, Js::AbstractArray, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)
    D = dim(model)
    invN = 1.0/nsites
    beta = 1.0/T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += model.cosines[model.spins[s]]
        M[2] += model.sines[model.spins[s]]
    end

    @inbounds for b in bonds(model)
        i, j = source(b), target(b)
        dir = bonddirection(b)
        dt = mod1(model.spins[j] - model.spins[i], model.Q)
        E -= model.cosines[dt] * Js[bondtype(b)]

        for d in 1:D
            U1[d] += model.cosines[dt] * dir[d]^2
            U2[d] += model.sines[dt] * dir[d]
        end
    end
    for d in 1:2
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    E /= nsites

    res = Measurement()
    for (i,c) in enumerate(["x","y"])
        res["Magnetization $c"] = M[i]
        res["|Magnetization $c|"] = abs(M[i])
        res["Magnetization $c^2"] = M[i]^2
        res["Magnetization $c^4"] = M[i]^4
        res["Helicity Modulus $c"] = U1[i]
    end
    M2 = sum(abs2, M)
    res["|Magnetization|"] = sqrt(M2)
    res["|Magnetization|^2"] = M2
    res["|Magnetization|^4"] = M2^2
    res["Energy"] = E
    res["Energy^2"] = E^2

    return res
end

function simple_estimator(model::XY, T::Real, Js::AbstractArray, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)
    D = dim(model)
    invN = 1.0/nsites
    beta = 1.0/T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += cospi(2model.spins[s])
        M[2] += sinpi(2model.spins[s])
    end

    @inbounds for b in bonds(model)
        i, j = source(b), target(b)
        dir = bonddirection(b)
        dt = mod(model.spins[j] - model.spins[i] + 2.0, 1.0)
        dt = ifelse(0.5 <= dt, dt - 1.0, dt)
        E -= cospi(2dt) * Js[bondtype(b)]

        for d in 1:D
            U1[d] += cospi(2dt) * dir[d]^2
            U2[d] += sinpi(2dt) * dir[d]
        end
    end
    for d in 1:2
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    E /= nsites

    res = Measurement()
    for (i,c) in enumerate(["x","y"])
        res["Magnetization $c"] = M[i]
        res["|Magnetization $c|"] = abs(M[i])
        res["Magnetization $c^2"] = M[i]^2
        res["Magnetization $c^4"] = M[i]^4
        res["Helicity Modulus $c"] = U1[i]
    end
    M2 = sum(abs2, M)
    res["|Magnetization|"] = sqrt(M2)
    res["|Magnetization|^2"] = M2
    res["|Magnetization|^4"] = M2^2
    res["Energy"] = E
    res["Energy^2"] = E^2

    return res
end

