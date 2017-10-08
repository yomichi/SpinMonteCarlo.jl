"""
    simple_estimate(model::Ising, T::Real, Js::AbstractArray)
    simple_estimate(model::Potts, T::Real, Js::AbstractArray)

return magnetization density `M` and total energy `E`.
"""
function simple_estimate(model::Ising, T::Real, Js::AbstractArray)
    nsites = numsites(model)
    nbonds = numbonds(model)

    M = mean(model.spins)
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(model, b), target(model, b)
        E += ifelse(model.spins[s1] == model.spins[s2], -1.0, 1.0) * Js[bondtype(model,b)]
    end
    E /= nsites
    return M,E
end

function simple_estimate(model::Potts, T::Real, Js::AbstractArray)
    nsites = numsites(model)
    nbonds = numbonds(model)
    invQ = 1.0/model.Q

    M = 0.0
    @inbounds for s in 1:nsites
        M += ifelse(model.spins[s]==1, 1.0-invQ, -invQ)
    end
    M /= nsites
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(model, b), target(model, b)
        E -= ifelse(model.spins[s1] == model.spins[s2], 1.0, 0.0) * Js[bondtype(model,b)]
    end
    E /= nsites
    return M,E
end


"""
    simple_estimate(model::Clock, T::Real, Js::AbstractArray)
    simple_estimate(model::XY, T::Real, Js::AbstractArray)

return magnetization density `M`, total energy `E`, and helicity modulus `U`.
"""
function simple_estimate(model::Clock, T::Real, Js::AbstractArray)
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

    @inbounds for b in 1:nbonds
        i, j = source(model, b), target(model,b)
        dir = bonddirection(model, b)
        dt = mod1(model.spins[j] - model.spins[i], model.Q)
        E -= model.cosines[dt] * Js[bondtype(model, b)]

        for d in 1:D
            U1[d] += model.cosines[dt] * dir[d]^2
            U2[d] += model.sines[dt] * dir[d]
        end
    end
    for d in 1:D
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    E /= nsites
    return M, E, U1
end

function simple_estimate(model::XY, T::Real, Js::AbstractArray)
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

    @inbounds for b in 1:nbonds
        i, j = source(model, b), target(model,b)
        dir = bonddirection(model, b)
        dt = mod(model.spins[j] - model.spins[i] + 2.0, 1.0)
        dt = ifelse(0.5 <= dt, dt - 1.0, dt)
        E -= cospi(2dt) * Js[bondtype(model,b)]

        for d in 1:D
            U1[d] += cospi(2dt) * dir[d]^2
            U2[d] += sinpi(2dt) * dir[d]
        end
    end
    for d in 1:D
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    E /= nsites
    return M, E, U1
end

